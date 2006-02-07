
#ifndef INCLUDED_MESSAGE_H
#define INCLUDED_MESSAGE_H

class Message
{
public:
  enum level { warning = 0, error=1 };
  
  void add(level l, std::string message)
  {
    if (l) m_errors.push_back(message);
    else m_warnings.push_back(message);
  }
  
  void display(std::ostream & os = std::cout)
  {
    for(unsigned int i=0; i<m_errors.size(); ++i)
      os << "ERROR:   " << m_errors[i] << std::endl;
    for(unsigned int i=0; i<m_warnings.size(); ++i)
      os << "WARNING: " << m_warnings[i] << std::endl;
  }
  
private:
  std::vector<std::string> m_errors;
  std::vector<std::string> m_warnings;
    
};

#endif
