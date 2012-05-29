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

#define FOR(i,n) for(int i=0, n_=(n); i<n_; i++)
#define MALLOC(type, size) (type*) malloc((size)*sizeof(type))
#define FREE(var)				\
  do {						\
    if(var) free(var);				\
    (var) = NULL;				\
  } while(0)

using namespace std;

typedef double Real;
typedef vector<Real> Vector;
typedef int Vertex;

/**
 * Some utility functions to work with vector
 */

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

inline bool operator == (const Vector &x, const Vector &y){
  FOR(i,x.size()) if(fabs(x[i]-y[i])>1e-10) return false;
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

// report error message and terminate the program 
inline void error_handle(const char* message){
  cout<<message<<endl;
  cout<<"Terminating the program"<<endl;
  exit(1);
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
class Timer{
private:
  HANDLE hProcess;
  ULONG64 start, end;
  double freq;
  
public:
  Timer(){
    DWORD processID = GetCurrentProcessId();
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
			                       PROCESS_VM_READ,
			                       FALSE, 
			                       processID );
             
    QueryProcessCycleTime(hProcess, &start);
	freq = 2.6e9;
  }
  
  double elapsed() {           
    QueryProcessCycleTime(hProcess, &end);
    
    ULONG64 diff = end - start;
    return double(diff)/freq;
  }
};

class Timer2{
private:
  HANDLE hProcess;
  ULONG64 start, end;
  double freq;
  
public:
  Timer2(){
    DWORD processID = GetCurrentProcessId();
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
			                       PROCESS_VM_READ,
			                       FALSE, 
			                       processID );
             
    QueryProcessCycleTime(hProcess, &start);
    
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
  
  double elapsed() { 
    DWORD processID = GetCurrentProcessId();
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
			                       PROCESS_VM_READ,
			                       FALSE, 
			                       processID );
             
    QueryProcessCycleTime(hProcess, &end);
    
    ULONG64 diff = end - start;
    return double(diff) / freq;
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

class Timer{
private:
  
  double start;
  
public:
  Timer(){
    struct timespec st;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &st);
    start = st.tv_sec + 1.0*st.tv_nsec/1e9;
  }
  
  double elapsed() { 
    struct timespec et;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &et);
    double end = et.tv_sec + 1.0*et.tv_nsec/1e9;
    return end - start;
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

