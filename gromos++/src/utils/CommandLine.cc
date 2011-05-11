#include "CommandLine.h"
#include "../../config.h"

#include <algorithm>

#if defined(HAVE_READLINE) && defined(HAVE_LIBREADLINE)
#include <cstdio>
#include <cstdlib>

#include <readline/readline.h>
#include <readline/history.h>
#endif

std::string utils::CommandLine::getLine(const std::string& prompt, std::ostream & out) {
  std::string line;
#if defined(HAVE_READLINE) && defined(HAVE_LIBREADLINE)
  char * buf;
  //rl_bind_key('\t',rl_abort);//disable auto-complete
  buf = readline(prompt.c_str());
  if (buf == NULL) return line;
  line = buf;
  if (buf[0]!=0) add_history(buf);
  free(buf);
#else
  if (!prompt.empty()) out << prompt;
  std::getline(std::cin, line);
#endif
  return line;
}

char to_lower(char c) { return std::tolower(c); }

bool utils::CommandLine::getYesNo(const std::string& prompt, std::ostream& out, const std::string& error) {
  while(true) {
    std::string a = getLine(prompt, out);
    std::transform(a.begin(), a.end(), a.begin(), to_lower);
    if (a[0] == 'y') return true;
    if (a[0] == 'n') return false;
    out << error;
  }
}
