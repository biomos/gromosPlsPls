#include <cassert>
#include <iomanip>
#include <algorithm>
#include "../args/Arguments.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "Hbond_calc_3c.h"
#include "Hbond_calc_2c.h"
#include "Hbond.h"

using utils::HB3c_calc;
using utils::HB3c;
using utils::Key2c;
using utils::Key3c;


void HB3c_calc::setval(const HB2c_calc& hb2c_calc, gcore::System& _sys, args::Arguments& _args) {
  sys = &_sys;
  args = &_args;
  initialise(hb2c_calc);
  set_reduce();
  pbc = args::BoundaryParser::boundary(_sys, _args);
  opents("Hbond_3c_time_index.out", "Hbond_3c_time_numHb.out");
}//end HB3c_calc::setval()

void HB3c_calc::clear() {
  --frames;
  hb3cc.clear();
  hb3cc_tmp.clear();
  ts.clear();
}//end HB3c_calc::clear()

void HB3c_calc::store_index(){
    native_key_storage.clear();
    for (HB3cContainer::const_iterator it = hb3cc_tmp.begin(); it != hb3cc_tmp.end() ; ++it) //go through native hbond map
    	native_key_storage.push_back(it->first);
}

void HB3c_calc::calc_native() {
  init_calc();

  for (size_t t = 0; t < native_key_storage.size(); ++t) {
    Key3c native_index = native_key_storage[t];

    int i,j,k,l;
    native_index.get_atom_num(i,j,l,k);
    std::vector<int> k_vec(1, k); //k must be a vector. in this case only with 1 element

    calc(i, j, k_vec);
  }
  ts.back().set_num(numHb);
}//end HB3c_calc::calc_native()

void HB3c_calc::calc(int i, int j, const std::vector<int>& k_atoms) {

  double angle1, dist1;
  gmath::Vec acceptor_j, bound_i, vec_ij;

  bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i), sys->box());
  acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j), sys->box());
  vec_ij = acceptor_j - donors.pos(i);

  if (neighbour(i, j) && distances(dist1, vec_ij) && angle(i, angle1, bound_i, vec_ij)) {

	std::vector<int>::const_iterator k_it = std::find(k_atoms.begin(), k_atoms.end(), j);
	if(k_it == k_atoms.end()) //j was not found in atoms_k
		k_it = k_atoms.begin();

    for(; k_it != k_atoms.end(); ++k_it){

      int k = *k_it;

      double angle2, dist2;
      gmath::Vec acceptor_k, vec_ik;

      acceptor_k = pbc->nearestImage(donors.pos(i), acceptors.pos(k), sys->box());
      vec_ik = acceptor_k - donors.pos(i);

		//       vvvv the inverse key could exist
      if (!hb3cc_tmp.count(Key3c(i,k,i,j)) && j != k && neighbour(i, k) && distances(dist2, vec_ik) && angle(i, angle2, bound_i, vec_ik)) {

        double angle_sum, dihedral;
        angle_sum = angle1 + angle2;

        if( anglesum(i, angle_sum, acceptor_j, acceptor_k) &&
            dihedrals(i, dihedral, bound_i, acceptor_j, acceptor_k)) {

            Key3c key(i,j,i,k); //d1,a1,d2,a2; d1==d2
            hb3cc_tmp[key].add(); //needed for reduce

            if(reduce){

                if(donors.mol(i) < 0){ //if atom i is from a solvent molecule
                    i = donors.atom(i) % donors.numSolventAtoms(); //get the solvent atom number (0,1,2,...)
                    int t;
                    for(t = 0; t < solv_donor.size()-1 && i != solv_donor[t]; ++t); //find the atom number of i: this is stored in solv_donor[t]. t = the position in donors
                    i = solv_donor.back() + t; //solv_donor.back stores where the first solvent starts in donors: add position of the solvent in donors to where the first solvent starts in donors
                }
                if(acceptors.mol(j) < 0){ //same for acceptors
                    j = acceptors.atom(j) % acceptors.numSolventAtoms();
                    int t;
                    for(t = 0; t < solv_acc.size()-1 && j != solv_acc[t]; ++t);
                    j = solv_acc.back() + t;
                }
                if(acceptors.mol(k) < 0){
                    k = acceptors.atom(k) % acceptors.numSolventAtoms();
                    int t;
                    for(t = 0; t < solv_acc.size()-1 && k != solv_acc[t]; ++t);
                    k = solv_acc.back() + t;
                }
                key.set_index(i,j,i,k);

                hb3cc[key].add(dist1, dist2, angle1, angle2, angle_sum, dihedral);
                ts.back().add_once(key);
            }
            else{
                hb3cc[key].add(dist1, dist2, angle1, angle2, angle_sum, dihedral);
                ts.back().add_all(key);
            }

            ++numHb;
          }
      }
    }
  }
}//end HB3c_calc::calc()


void HB3c_calc::calc_hb(CubeSystem<int>& cubes_donors, CubeSystem<int>& cubes_acceptors){
    init_calc();
//#################

    if(cubes_donors.size() != cubes_acceptors.size())
        throw gromos::Exception("hbond","Donor & Acceptor CubeSystems do not have the same size. Please treat them equally!");
    cubes_donors.delete_atoms();
	cubes_acceptors.delete_atoms();
	//donX && accY:
    if(!(donX.empty() || accY.empty())){

		//assign donX atoms to cubes
		for (int i = 0; i < donX.size(); ++i)
		    cubes_donors.assign_atom(donX[i], donors.pos(donX[i]));

		//assign accY atoms to cubes
		for (int i = 0; i < accY.size(); ++i)
		    cubes_acceptors.assign_atom(accY[i], acceptors.pos(accY[i]));

		if(cubes_donors.size() != cubes_acceptors.size())
		    throw gromos::Exception("hbond","Donor & Acceptor CubeSystems do not have the same size. Please treat them equally!");

		for(int c=0; c<cubes_donors.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes_donors and cubes_acceptors

		    int don_atoms_i=cubes_donors.cube_i_atomlist(c).size(), //donor atoms my cube
		        don_atoms_j=cubes_donors.cube_j_atomlist(c).size(), //donor atoms neighbour cube
		        acc_atoms_i=cubes_acceptors.cube_i_atomlist(c).size(), //acceptor atoms my cube
		        acc_atoms_j=cubes_acceptors.cube_j_atomlist(c).size(); //acceptor atoms neighbour cube

			std::vector<int> k_list;
			if(don_atoms_i && acc_atoms_j || don_atoms_j && acc_atoms_i)
                k_list = cubes_acceptors.neighbour_atomlist(c); //acceptor atoms k list

		    if(don_atoms_i && acc_atoms_j) {
		    //no skiP: we do not compare the same atoms in i and j!
		        for(int di=0; di < don_atoms_i; ++di ){ //go through all donor atoms in first cube...
		            for(int aj=0; aj < acc_atoms_j; ++aj){ //..and check against all acceptor atoms in the neighbour cube
		                calc(cubes_donors.cube_i_atomlist(c).at(di), cubes_acceptors.cube_j_atomlist(c).at(aj), k_list);
		            }
		        }
		    }
		    if(don_atoms_j && acc_atoms_i && !cubes_donors.same_cube(c)){
		        for(int dj=0; dj < don_atoms_j; ++dj ){ //then do the same for the donor atoms in the second cube. this must be done since the cubes do not have redundancy. every cube pair only exists once. e.g. only once instance of 0-2 and 2-0
		            for(int ai=0; ai < acc_atoms_i; ++ai){
		                calc(cubes_donors.cube_j_atomlist(c).at(dj), cubes_acceptors.cube_i_atomlist(c).at(ai), k_list);
		            }
		        }
		    }
		}
    }
//#################
    //donZ & accX+accY:
    if(!(donZ.empty() || accX.empty() && accY.empty() )){ //if there is something in accB

        cubes_donors.delete_atoms(); //delete atoms from cubes
        //we can keep accY atoms and add accX

        //assign donZ
        for (int i = 0; i < donZ.size(); ++i)
            cubes_donors.assign_atom(donZ[i], donors.pos(donZ[i]));

        //assign accX atoms to cubes
        for (int i = 0; i < accX.size(); ++i)
            cubes_acceptors.assign_atom(accX[i], acceptors.pos(accX[i]));

        for(int c=0; c<cubes_donors.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes_donors and cubes_acceptors

            int don_atoms_i=cubes_donors.cube_i_atomlist(c).size(), //donor atoms my cube
		        don_atoms_j=cubes_donors.cube_j_atomlist(c).size(), //donor atoms neighbour cube
		        acc_atoms_i=cubes_acceptors.cube_i_atomlist(c).size(), //acceptor atoms my cube
		        acc_atoms_j=cubes_acceptors.cube_j_atomlist(c).size(); //acceptor atoms neighbour cube

			std::vector<int> k_list;
			if(don_atoms_i && acc_atoms_j || don_atoms_j && acc_atoms_i)
                k_list = cubes_acceptors.neighbour_atomlist(c); //acceptor atoms k list

			if(don_atoms_i && acc_atoms_j) {
            //no skiP: we do not compare the same atoms in i and j!
                for(int di=0; di < don_atoms_i; ++di ){ //go through all donor atoms in first cube...
                    for(int aj=0; aj < acc_atoms_j; ++aj){ //..and check against all acceptor atoms in the neighbour cube
                        calc(cubes_donors.cube_i_atomlist(c).at(di), cubes_acceptors.cube_j_atomlist(c).at(aj), k_list);
                    }
                }
            }
            if(don_atoms_j && acc_atoms_i && !cubes_donors.same_cube(c)){
                for(int dj=0; dj < don_atoms_j; ++dj ){ //then do the same for the donor atoms in the second cube. this must be done since the cubes do not have redundancy. every cube pair only exists once. e.g. only once instance of 0-2 and 2-0
                    for(int ai=0; ai < acc_atoms_i; ++ai){
                        calc(cubes_donors.cube_j_atomlist(c).at(dj), cubes_acceptors.cube_i_atomlist(c).at(ai), k_list);
                    }
                }
            }
        }
    }
//#########
//donX & donY + accX & accZ
//this will always be called

	cubes_donors.delete_atoms(); //delete atoms from cubes
    cubes_acceptors.delete_atoms();

    //assign donX
    for (int i = 0; i < donX.size(); ++i)
        cubes_donors.assign_atom(donX[i], donors.pos(donX[i]));
    //assign donY
    for (int i = 0; i < donY.size(); ++i)
        cubes_donors.assign_atom(donY[i], donors.pos(donY[i]));

    //assign accX atoms to cubes
    for (int i = 0; i < accX.size(); ++i)
        cubes_acceptors.assign_atom(accX[i], acceptors.pos(accX[i]));
    //assign accZ atoms to cubes
    for (int i = 0; i < accZ.size(); ++i)
        cubes_acceptors.assign_atom(accZ[i], acceptors.pos(accZ[i]));

    for(int c=0; c<cubes_donors.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes_donors and cubes_acceptors

        int don_atoms_i=cubes_donors.cube_i_atomlist(c).size(), //donor atoms my cube
            don_atoms_j=cubes_donors.cube_j_atomlist(c).size(), //donor atoms neighbour cube
            acc_atoms_i=cubes_acceptors.cube_i_atomlist(c).size(), //acceptor atoms my cube
            acc_atoms_j=cubes_acceptors.cube_j_atomlist(c).size(); //acceptor atoms neighbour cube

        std::vector<int> k_list;
        if(don_atoms_i && acc_atoms_j || don_atoms_j && acc_atoms_i)
            k_list = cubes_acceptors.neighbour_atomlist(c); //acceptor atoms k list

        if(don_atoms_i && acc_atoms_j) {
        //no skiP: we do not compare the same atoms in i and j!
            for(int di=0; di < don_atoms_i; ++di ){ //go through all donor atoms in first cube...
                for(int aj=0; aj < acc_atoms_j; ++aj){ //..and check against all acceptor atoms in the neighbour cube
                    calc(cubes_donors.cube_i_atomlist(c).at(di), cubes_acceptors.cube_j_atomlist(c).at(aj), k_list);
                }
            }
        }
        if(don_atoms_j && acc_atoms_i && !cubes_donors.same_cube(c)){
            for(int dj=0; dj < don_atoms_j; ++dj ){
                for(int ai=0; ai < acc_atoms_i; ++ai){
                    calc(cubes_donors.cube_j_atomlist(c).at(dj), cubes_acceptors.cube_i_atomlist(c).at(ai), k_list);
                }
            }
        }
    }
    ts.back().set_num(numHb);
}

void HB3c_calc::calc_vac(){
    init_calc();

	  // donX + accY
    if(!(donX.empty() || accY.empty()))
        for (int i = 0; i < donX.size(); ++i)
            for (int j = 0; j < accY.size(); ++j)
                calc(donX[i], accY[j], accY);

	//donZ + accX +accY
	if(!(donZ.empty() || (accX.empty() && accY.empty()) ))
        for (int i = 0; i < donZ.size(); ++i){
            for (int j = 0; j < accX.size(); ++j)
                calc(donZ[i], accX[j], accX);
			for (int j = 0; j < accY.size(); ++j)
                calc(donZ[i], accY[j], accY);
		}

	//donX + accX +accZ
    if(!(donX.empty() || (accX.empty() && accZ.empty()) ))
        for (int i = 0; i < donX.size(); ++i){
            for (int j = 0; j < accX.size(); ++j)
                calc(donX[i], accX[j], accX);
            for (int j = 0; j < accZ.size(); ++j)
                calc(donX[i], accZ[j], accZ);
        }
	//donY + accX +accZ
	if(!(donY.empty() || (accX.empty() && accZ.empty()) ))
        for (int i = 0; i < donY.size(); ++i){
            for (int j = 0; j < accX.size(); ++j)
                calc(donY[i], accX[j], accX);
            for (int j = 0; j < accZ.size(); ++j)
                calc(donY[i], accZ[j], accZ);
        }
    ts.back().set_num(numHb);
}

void HB3c_calc::printstatistics(bool sort_occ, double higher){

  cout << endl << "# Three-centered hydrogen bonds:" << endl;

    print_header();

    unsigned int i = 1;
    for(HB3cContainer::iterator it = hb3cc.begin(); it != hb3cc.end(); ++it, ++i){
        it->second.set_id(i);
        if(it->second.num()/(double)frames * 100 >= higher)
            print(it->first);
	}

    std::sort(ts.begin(), ts.end(), CompTime<Key3c>);

    //write timeseries - numhb and   timeseries- hbindex
    for(TimeseriesContainer::const_iterator it = ts.begin(); it != ts.end(); ++it){

        timeseriesHBtot << setw(15) << it->time()
                        << setw(10) << it->num() << endl;

        for(std::vector<Key3c>::const_iterator it_key = it->keys().begin(); it_key != it->keys().end(); ++it_key){
            timeseriesHB << setw(15) << it->time()
                         << setw(10) << hb3cc[*it_key].id() << endl;
        }
    }

    if(sort_occ){
        HB3cContainerIteratorList hb_vec; //vector of iterators to hb2cc map

	    for(HB3cContainer::iterator it = hb3cc.begin(); it != hb3cc.end(); ++it){
    	    hb_vec.push_back(it); //populate vector with iterators to map
		}

        std::sort(hb_vec.begin(), hb_vec.end(), sort_rev_by_occ<Key3c, HB3c>);

        cout << endl << "# SORTED three-centered hydrogen bonds:" << endl;
        print_header();

		for(HB3cContainerIteratorList::const_iterator it = hb_vec.begin(); it!= hb_vec.end(); ++it){
            if((**it).second.num()/(double)frames * 100 < higher) // as soon as the occurence drops below higher value: stop
                break;
        	print((**it).first);
        }
    }

}//end HB3c_calc::printstatistics()

void HB3c_calc::print_header() const{
    cout  << right
          << "#"
          << setw(8) << "HB"
          << setw(16) << "Donor"
          << setw(20) << "Acceptor"
          << setw(18) << "D -"
          << setw(15) << "H ..."
          << setw(11) << "A"
          << setw(13) << "DIST"
          << setw(8) << "ANGLE"
          << setw(9) << "SUM"
          << setw(8) << "DIHED."
          << setw(9) << "OCCUR"
          << setw(11) << "%"
          << endl;

}

void HB3c_calc::print(const Key3c& key){

    int i_d, i_a1, i_a2,i_d2;
    key.get_atom_num(i_d,i_a1,i_d2,i_a2);
    const HB3c& hb3cprint = hb3cc[key]; //= HB3c
    int occur = hb3cprint.num(); //how many hbonds of this index

    cout << setw(8) << hb3cc[key].id();

    if (donors.mol(i_d) < 0) cout << setw(8) << " ";
    else cout << setw(8) << donors.mol(i_d) + 1;
    cout << setw(5) << donors.resnum(i_d) + 1
          << setw(6) << donors.resname(i_d)
          << setw(2) << "-";
    if (acceptors.mol(i_a1) < 0) cout << setw(4) << " ";
    else cout << setw(4) << acceptors.mol(i_a1) + 1;
    cout << setw(5) << acceptors.resnum(i_a1) + 1
          << setw(6) << acceptors.resname(i_a1)
          << setw(11) << bound.atom(i_d) + 1
          << setw(6) << bound.name(i_d)
          << setw(2) << "-"
          << setw(6) << donors.atom(i_d) + 1
          << setw(6) << donors.name(i_d)
          << setw(2) << "-"
          << setw(6) << acceptors.atom(i_a1) + 1
          << setw(6) << acceptors.name(i_a1)
          << setprecision(3) << setw(13) << hb3cprint.meandist(0)
          << setw(8) << hb3cprint.meanangle(0)
          << std::fixed
          << setprecision(3) << setw(9) << hb3cprint.meanangle_sum()
          << setw(8) << hb3cprint.meandihedral() << " "
          << setprecision(0) << setw(8) << hb3cprint.num() << " "
          << setprecision(2) << setw(10) << ((hb3cprint.num() / (double) frames)*100)
          << endl;

    // and the second line
    cout << setw(27) << " "
          << setw(2) << "-";
    if (acceptors.mol(i_a2) < 0) cout << setw(4) << " ";
    else cout << setw(4) << acceptors.mol(i_a2) + 1;
    cout << setw(5) << acceptors.resnum(i_a2) + 1
          << setw(6) << acceptors.resname(i_a2)
          << setw(27) << " "
          << setw(6) << "-"
          << setw(6) << acceptors.atom(i_a2) + 1
          << setw(6) << acceptors.name(i_a2)
          << setprecision(3) << setw(13) << hb3cprint.meandist(1)
          << setw(8) << hb3cprint.meanangle(1)
          << endl;
}

//merge hbond objects
void HB3c_calc::merge(utils::HB3c_calc& input){
    #ifdef OMP
    #pragma omp critical
    #endif
    {   //merge maps:
        for(HB3cContainer::const_iterator it=input.hb3cc.begin(); it != input.hb3cc.end(); ++it){ //go through input map
                hb3cc[it->first].add(it->second); //add the HB3c entries
        }

        //merge other things:
        frames += input.frames; //add frames so we get the total number of frames
        ts.insert(ts.end(), input.ts.begin(), input.ts.end());
      }
}

