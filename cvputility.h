#ifndef __CVPUTILITY_H__
#define __CVPUTILITY_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <cstdarg>
#include <map>
#include <typeinfo>
#include <cctype>
#include <algorithm>
#include <utility>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Sparse>

#define FOR(i,n) for(int i=0, n_=(n); i<n_; i++)
#define MALLOC(type, size) (type*) malloc((size)*sizeof(type))
#define ITER(vec, it) for(Vector::InnerIterator it(vec); it; ++it)
#define FREE(var)				\
  do {						\
    if(var) free(var);				\
    (var) = NULL;				\
  } while(0)

using namespace std;
using namespace Eigen;

typedef double Real;
typedef SparseVector<Real> Vector;
//typedef vector<Real> Vector;
typedef int Vertex;

// report error message and terminate the program 
inline void error_handle(const string &message){
  cout<<message<<endl;
  cout<<"Terminating the program"<<endl;
  exit(1);
}

class TableReport {
 private:
  string row_format, header_format, line;

 public:
  TableReport(const string &format) {
    string row_token, header_token;
    string specifiers = "cdieEfgGosuxXpn";
    bool row_token_flag = false, header_token_flag = false;

    line = row_format = header_format = row_token = header_token = "";

    FOR(i, format.length()){
      if(format[i] == '%') 
	row_token_flag = header_token_flag = true;
      if(row_token_flag){
	row_token += format[i];
	if(row_token_flag && specifiers.find(format[i]) != string::npos){
	  row_format += row_token;
	  row_token = "";
	  row_token_flag = false;
	}
      }
      
      if(header_token_flag){
	if(format[i] == '.' || isalpha(format[i])){
	  header_format += header_token + "s";
	  header_token = "";
	  header_token_flag = false;
	}
	else header_token += format[i];
      }
    }
  }

  fstream& print_header(fstream &f, ...){
    char header[300];
    va_list args;
    va_start(args, f);
    vsprintf(header, header_format.c_str(), args);

    line = "";
    FOR(i, strlen(header)) line += "-"; 

    f<<line<<endl<<header<<endl<<line<<endl;

    va_end(args);
    return f;
  }

  fstream& print_row(fstream &f, ...){
    char row[300];
    va_list args;
    va_start(args, f);
    vsprintf(row, row_format.c_str(), args);
    f<<row<<endl;
    va_end(args);
    return f;
  }

  fstream& print_line(fstream &f){
    f<<line<<endl;
    return f;
  }
};

class Timer{
 private:
  vector<double> records;

 protected:
  virtual double get_clock() const{
    return double(clock())/double(CLOCKS_PER_SEC);
  }

 public:
  void record(){
    records.push_back(get_clock());
  }

  double elapsed(int i=-1, int j=-2){
    int n = records.size();
    if(n<2) return 0.0;
    if(i<0) i+=n;
    if(j<0) j+=n;
    return fabs(records[i]-records[j]);
  }
};

class SettingMapper{
 private:
  map<string, int> int_params;
  map<string, double> real_params;
  map<string, string> str_params;

  string trim(const string &key){
    string punc = "+-.\\/";
    char mark = '$';
    string skey(key);
    int n = skey.length();

    FOR(i, n) 
      if(!isalnum(skey[i]) && punc.find(skey[i])==string::npos) 
	skey[i] = mark; 
      else break;
    FOR(i, n) 
      if(!isalnum(skey[n-1-i]) && punc.find(skey[n-1-i])==string::npos) 
	skey[n-1-i] = mark; 
      else break;
    int count_marks = count(skey.begin(), skey.end(), mark);
    remove(skey.begin(), skey.end(), mark);
    skey.erase(skey.length() - count_marks);
    return skey;
  }

  string tolower_str(const string &key){
    string skey(key);
    FOR(i, skey.length()) skey[i] = tolower(skey[i]);
    return skey;
  }

 public:
  
  void seti(const string &key, int value){ 
    int_params[tolower_str(trim(key))] = value;
  }

  void setr(const string &key, double value){ 
    real_params[tolower_str(trim(key))] = value;
  }

  void sets(const string &key, const string &value){ 
    str_params[tolower_str(trim(key))] = value;
  }

  int& geti(const string &key){
    string skey = tolower_str(trim(key));
    if(int_params.find(skey) == int_params.end())
      error_handle("Setting "+key+" of type int has not been set");
    return int_params[skey];
  }  

  double& getr(const string &key){
    string skey = tolower_str(trim(key));
    if(real_params.find(skey) == real_params.end())
      error_handle("Setting "+key+" of type real has not been set");
    return real_params[skey];
  }  

  string& gets(const string &key){
    string skey = tolower_str(trim(key));
    cout<<skey<<endl;
    if(str_params.find(skey) == str_params.end())
      error_handle("Setting "+key+" of type string has not been set");
    cout<<str_params[skey]<<endl;
    return str_params[skey];
  }

  void set_by_str(const string &key, const string &value_str){
    char *pend = NULL;
    string svalue = trim(value_str);
    string skey = tolower_str(trim(key));

    long int int_value = strtol(svalue.c_str(), &pend, 10);
    if(strlen(pend) == 0){
      cout<<"Setting "<<skey<<"="<<int_value<<endl;
      int_params[skey] = int(int_value);
      return;
    }

    double real_value = strtod(svalue.c_str(), &pend);
    if(strlen(pend) == 0){
      cout<<"Setting "<<skey<<"="<<real_value<<endl;
      real_params[skey] = real_value;
      return;
    }
    
    if(svalue == "yes" || svalue == "no")
      int_params[skey] = int(svalue == "yes");
    else{
      cout<<"Setting "<<skey<<"="<<svalue<<endl;
      str_params[skey] = svalue;
    }
  }

  void read(fstream &f){
    string line, key, value;
    while(!f.eof()){
      getline(f, line);
      size_t delim = line.find('=');
      key = line.substr(0, delim);
      value = line.substr(delim+1);
      set_by_str(key, value);
    }
  }
};


/**
 * Some utility functions to work with vector
 */

/*
inline Vector operator + (const Vector &x, const Vector &y){
  Vector z(x);
  FOR(i, x.size()) z[i] += y[i];
  return z;
}

inline Vector operator - (const Vector &x, const Vector &y){
  Vector z(x);
  FOR(i, x.size()) z[i] -= y[i];
  return z;
}

inline Vector & operator += (Vector &x, const Vector &y){
  FOR(i, x.size()) x[i] += y[i];
  return x;
}

inline Vector & operator -= (Vector &x, const Vector &y){
  FOR(i, x.size()) x[i] -= y[i];
  return x;
}

inline Vector & operator *= (Vector &x, Real alpha){
  FOR(i, x.size()) x[i] *= alpha;
  return x;
}

inline Real operator * (const Vector &x, const Vector &y){
  Real sum = 0.0;
  FOR(i, x.size()) sum += x[i]*y[i];
  return sum;
}

inline Vector operator *(Real alpha, const Vector &x){
  Vector z(x);
  z*=alpha;
  return z;
}

*/

inline bool operator == (const Vector &x, const Vector &y){
  Vector::InnerIterator itx(x), ity(y);
  for(; itx && ity; ++itx, ++ity){
    if(itx.index() != ity.index()) return false;
    if(fabs(itx.value()-ity.value())>1e-10) return false;
  }
  if(itx) return false;
  if(ity) return false;
  return true;
}


template <typename T> inline bool updatemin(T &a, T b){
  if(a>b) return a=b, true;
  return false;
}

template <typename T> inline bool updatemax(T &a, T b){
  if(a<b) return a=b, true;
  return false;
}

template<typename T>
bool read_param(string fininame, const string &alias, T &t){
  fstream f(fininame.c_str(), fstream::in);
  char line[305], *str;

  if(!f) error_handle("No CVP.ini file found.");
  while(!f.eof()){
    f.getline(line, 300);
    str = strrchr(line,'#'); // comments are marked by '#' character
    if(str != NULL) continue;
    if(strstr(line, alias.c_str()) == NULL) continue;
    str = strrchr(line,'='); // value delim by the '=' character
    if(str == NULL) continue; // no delim found, read next line

    stringstream ss(str+1, stringstream::in);
    ss>>t;
    f.close();
    return true;
  }

  f.close();
  return false;
}

/**
 * ========= Timing and Memory Reporting =============================
 *
 * Memory and CPU usage information is dependent on archi and platform
 * Therefore, two different timer and memory report are implemented
 * one for Windows and one for Linux
 *
 * Windows:
 *     Windows timer makes use of QuerryProcessCycleTime() (kernel32.lib) to 
 *     get the CPU cycles and CallNtPowerInformation() (PowrProf.lib) to get 
 *     the CPU frequency
 *     Windows memory queries from GetProcessMemoryInfo() (psapi.lib)
 *    
 *     requires: psapi.lib, PowrProf.lib
 *
 * Linux:
 *     Linux timer use clock_gettime() from /sys/times.h
 *     Linux memory report is read from /proc/self/status
 *     
 * ====================================================================
 */

#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#include <Powrprof.h>

//Declaration of CPUcycletime: 
class CPUTimer : public Timer{
 private:
  HANDLE hProcess;
  ULONG64 start, end;

 protected:
  double freq;

  virtual double get_clock(){
    ULONG64 end;
    QueryProcessCycleTime(hProcess, &end);
    return double(end-start)/freq;
  }
  
public:
  CPUTimer(double freg_ = 2.6e9) : Timer(), freq(freq_) {
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
			     PROCESS_VM_READ,
			     FALSE, 
			     GetCurrentProcessId() );
    QueryProcessCycleTime(hProcess, &start);
    freq = 2.6e9;
  }
  
};

class CPUTimer2 : public CPUTimer{
public:
  CPUTimer2() : CPUTimer() {
    typedef struct _PROCESSOR_POWER_INFORMATION {
      ULONG  Number;
      ULONG  MaxMhz;
      ULONG  CurrentMhz;
      ULONG  MhzLimit;
      ULONG  MaxIdleState;
      ULONG  CurrentIdleState;
    } PROCESSOR_POWER_INFORMATION, *PPROCESSOR_POWER_INFORMATION;

    SYSTEM_INFO si = {0};
    ::GetSystemInfo(&si);
    int size = si.dwNumberOfProcessors * sizeof(PROCESSOR_POWER_INFORMATION);
    LPBYTE pBuffer = new BYTE[size]; 
    NTSTATUS status = ::CallNtPowerInformation(ProcessorInformation, NULL, 0, pBuffer, size);
    PPROCESSOR_POWER_INFORMATION ppi = (PPROCESSOR_POWER_INFORMATION)pBuffer;
    freq = (ppi->CurrentMhz)*1e6;
  }  
};

inline void memory_report(fstream &f){
  DWORD processID = GetCurrentProcessId();
  HANDLE hProcess;
  PROCESS_MEMORY_COUNTERS pmc;
  // Print the process identifier.
  f<<endl<<"Memory report for Process ID "<<processID<<endl;

  // Print information about the memory usage of the process.
  hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
			   PROCESS_VM_READ,
			   FALSE, 
			   processID );
  if (NULL == hProcess) return;

  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
    //f<<"\tPageFaultCount: "<<pmc.PageFaultCount<<endl;
    f<<"PEAK MEMORY CONSUMPTION: "<<pmc.PeakWorkingSetSize<<endl;
    //f<<"\tYour app's CURRENT MEMORY CONSUMPTION: "<<pmc.WorkingSetSize<<endl;
    //f<<"\tQuotaPeakPagedPoolUsage: "<<pmc.QuotaPeakPagedPoolUsage<<endl;
    //f<<"\tQuotaPagedPoolUsage: "<<pmc.QuotaPagedPoolUsage<<endl;
    //f<<"\tQuotaPeakNonPagedPoolUsage: "<<pmc.QuotaPeakNonPagedPoolUsage<<endl;
    //f<<"\tQuotaNonPagedPoolUsage: "<<pmc.QuotaNonPagedPoolUsage<<endl;
    //f<<"\tPagefileUsage: "<<pmc.PagefileUsage<<endl; 
    //f<<"\tPeakPagefileUsage: "<<pmc.PeakPagefileUsage<<endl;
  }
  CloseHandle( hProcess );
}

#elif __linux__
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

class CPUTimer : public Timer{
 protected:
  virtual double get_clock() { 
    struct timespec et;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &et);
    return double(et.tv_sec) + 1.0e-9*double(et.tv_nsec);
  }
};

inline void memory_report(fstream &f){
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      stringstream ss(line+7, stringstream::in);
      ss>>result;
      f<<"Peak Memory Usage: "<< result << "kB" << endl;
      break;
    }
  }
  fclose(file);
}
#endif

#endif

