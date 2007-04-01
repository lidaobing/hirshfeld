#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "base_directory.hpp"
#include <vector>
#include "dirname.h"

using namespace std;

class BaseDirectory {
public:
  static const BaseDirectory* get_instance() {
    if(not m_instance) {
      m_instance = new BaseDirectory();
    }
    return m_instance;
  }

  const string& home() const {return m_home;}
  const string& xdg_data_home() const {return m_xdg_data_home;}
  const vector<string>& xdg_data_dirs() const {return m_xdg_data_dirs;}
private:
  BaseDirectory() {
    const char* tmp;

    tmp = getenv("HOME");
    if(tmp) {
      m_home = tmp;
    } else {
      m_home = "/";
    }

    tmp = getenv("XDG_DATA_DIRS");
    if(tmp) {
      m_xdg_data_home = tmp;
    } else {
      m_xdg_data_home = m_home + "/.local/share";
    }

    tmp = getenv("XDG_DATA_DIRS");
    if(not tmp) tmp = "/usr/local/share:/usr/share";

    string tmp2 = tmp;
    while(not tmp2.empty() and tmp2[0] == ':') {
      tmp2 = tmp2.substr(1);
    }

    while(not tmp2.empty() and tmp2[tmp2.size()-1] == ':') {
      tmp2 = tmp2.substr(0, tmp2.size()-1);
    }

    m_xdg_data_dirs.push_back(m_xdg_data_home);
    size_t pos = 0;
    while(pos != string::npos) {
      size_t nextpos = tmp2.find(':', pos);
      if(nextpos == string::npos) {
        m_xdg_data_dirs.push_back(tmp2.substr(pos));
        pos = nextpos;
      } else {
        m_xdg_data_dirs.push_back(tmp2.substr(pos, nextpos-pos));
        pos = nextpos+1;
      }
    }
  }
  
  string m_home;
  string m_xdg_data_home;
  vector<string> m_xdg_data_dirs;
  static BaseDirectory* m_instance;
};
BaseDirectory* BaseDirectory::m_instance = 0;

static string home;
static string xdg_data_home;
static vector<string> xdg_data_dirs;

std::string
load_first_data(const std::string& resource)
{
  const BaseDirectory& bd = *(BaseDirectory::get_instance());
  const vector<string>& xdg_data_dirs = bd.xdg_data_dirs();

  for(int i = 0; i < int(xdg_data_dirs.size()); ++i) {
    string res = xdg_data_dirs[i] + "/" + resource;
    if(access(res.c_str(), R_OK) == 0) return res;
  }
  return "";
}
