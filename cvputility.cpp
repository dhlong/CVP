#include "cvputility.h"

void error_handle(const string &message){
  cout<<message<<endl;
  cout<<"Terminating the program"<<endl;
  exit(1);
}

TableReport::TableReport(const string &format) {
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

ostream& TableReport::print_header(ostream *f, ...){
  char header[300];
  va_list args;
  va_start(args, f);
  vsprintf(header, header_format.c_str(), args);
  
  line = "";
  FOR(i, strlen(header)) line += "-"; 
  
  (*f)<<line<<endl<<header<<endl<<line<<endl;
  
  va_end(args);
  return *f;
}

ostream& TableReport::print_row(ostream *f, ...){
  char row[300];
  va_list args;
  va_start(args, f);
  vsprintf(row, row_format.c_str(), args);
  (*f)<<row<<endl;
  va_end(args);
  return *f;
}

ostream& TableReport::print_line(ostream &f){
  f<<line<<endl;
  return f;
}

double Timer::get_clock() const{
  return double(clock())/double(CLOCKS_PER_SEC);
}

void Timer::record(){
  records.push_back(this->get_clock());
}

double Timer::elapsed(int i, int j){
  int n = records.size();
  if(n<2) return 0.0;
  if(i<0) i+=n;
  if(j<0) j+=n;
  return fabs(records[i]-records[j]);
}


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

void SettingMapper::setb(const string &key, bool value){ 
  bool_params[tolower_str(trim(key))] = value;
}

void SettingMapper::seti(const string &key, int value){ 
  int_params[tolower_str(trim(key))] = value;
}

void SettingMapper::setr(const string &key, double value){ 
  real_params[tolower_str(trim(key))] = value;
}

void SettingMapper::sets(const string &key, const string &value){ 
  str_params[tolower_str(trim(key))] = value;
}

bool& SettingMapper::getb(const string &key){
  string skey = tolower_str(trim(key));
  if(bool_params.find(skey) == bool_params.end())
    error_handle("Setting '"+key+"' of type bool has not been set.");
  return bool_params[skey];
}  

int& SettingMapper::geti(const string &key){
  string skey = tolower_str(trim(key));
  if(int_params.find(skey) == int_params.end())
    error_handle("Setting '"+key+"' of type int has not been set.");
  return int_params[skey];
}  

double& SettingMapper::getr(const string &key){
  string skey = tolower_str(trim(key));
  if(real_params.find(skey) == real_params.end())
    error_handle("Setting '"+key+"' of type real has not been set.");
  return real_params[skey];
}  

string& SettingMapper::gets(const string &key){
  string skey = tolower_str(trim(key));
  cout<<skey<<endl;
  if(str_params.find(skey) == str_params.end())
    error_handle("Setting '"+key+"' of type string has not been set.");
  cout<<str_params[skey]<<endl;
  return str_params[skey];
}

void SettingMapper::set_by_str(const string &key, const string &value_str){
  char *pend = NULL;
  string svalue = trim(value_str);
  string skey = tolower_str(trim(key));

  //cout<<"           Setting "<<key<<" = "<<value_str<<endl;
  
  long int int_value = strtol(svalue.c_str(), &pend, 10);
  if(strlen(pend) == 0){
    int_params[skey] = int(int_value);
    return;
  }
  
  double real_value = strtod(svalue.c_str(), &pend);
  if(strlen(pend) == 0){
    real_params[skey] = real_value;
    return;
  }
  
  if(svalue == "yes" || svalue == "no")
    bool_params[skey] = (svalue == "yes");
  else
    str_params[skey] = svalue;
}

void SettingMapper::read(fstream &f){
  string line, key, value;
  cout<<"Reading settings"<<endl;
  while(!f.eof()){
    getline(f, line);
    size_t delim = line.find('=');
    key = line.substr(0, delim);
    value = line.substr(delim+1);
    // cout<<"********** Setting "<<key <<" = "<<value<<endl;
    set_by_str(key, value);
  }
}

void SettingMapper::report(fstream &f){
#define MAPITER(type, var, it)			\
  for(map<string, type>::iterator it = (var).begin(); it != (var).end(); ++it)

  MAPITER(bool, bool_params, it){
    f<<"\t"<<setw(50)<<left<<(*it).first<<" = ";
    f<< (((*it).second) ? "YES" : "NO") <<endl;
  }

  MAPITER(int, int_params, it)
    f<<"\t"<<setw(50)<<left<<(*it).first<<" = "<<(*it).second<<endl;

  MAPITER(double, real_params, it)
    f<<"\t"<<setw(50)<<left<<(*it).first<<" = "<<(*it).second<<endl;

  MAPITER(string, str_params, it)
    f<<"\t"<<setw(50)<<left<<(*it).first<<" = "<<(*it).second<<endl;
}

/*
bool operator == (const Vector &x, const Vector &y){
  Vector::InnerIterator itx(x), ity(y);
  for(; itx && ity; ++itx, ++ity){
    if(itx.index() != ity.index()) return false;
    if(fabs(itx.value()-ity.value())>1e-10) return false;
  }
  if(itx) return false;
  if(ity) return false;
  return true;
}
*/

#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#include <Powrprof.h>

//Declaration of CPUcycletime: 
double CPUTimer::get_clock() const{
  HANDLE hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
				  PROCESS_VM_READ,
				  FALSE, 
				  GetCurrentProcessId() );
  ULONG64 end;
  double freq = 2.6e9;

  QueryProcessCycleTime( hProcess,&end );
  CloseHandle(hProcess);
  return double(end)/freq;
}

/*  
CPUTimer2:CPUTimer2() : CPUTimer() {
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
*/ 

//
//f<<"\tPageFaultCount: "<<pmc.PageFaultCount<<endl;
//f<<"\tYour app's PEAK MEMORY CONSUMPTION: "<<pmc.PeakWorkingSetSize<<endl;
//f<<"\tYour app's CURRENT MEMORY CONSUMPTION: "<<pmc.WorkingSetSize<<endl;
//f<<"\tQuotaPeakPagedPoolUsage: "<<pmc.QuotaPeakPagedPoolUsage<<endl;
//f<<"\tQuotaPagedPoolUsage: "<<pmc.QuotaPagedPoolUsage<<endl;
//f<<"\tQuotaPeakNonPagedPoolUsage: "<<pmc.QuotaPeakNonPagedPoolUsage<<endl;
//f<<"\tQuotaNonPagedPoolUsage: "<<pmc.QuotaNonPagedPoolUsage<<endl;
//f<<"\tPagefileUsage: "<<pmc.PagefileUsage<<endl; 
//f<<"\tPeakPagefileUsage: "<<pmc.PeakPagefileUsage<<endl;
//
double memory_usage(int peak){
  HANDLE hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
				  PROCESS_VM_READ,
				  FALSE, 
				  GetCurrentProcessId() );
  PROCESS_MEMORY_COUNTERS pmc;

  if (NULL == hProcess) return -1.0; // error

  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
    if(peak == 1) return double(pmc.PeakWorkingSetSize);
    else return double(pmc.WorkingSetSize);
  else 
    return -1.0;

  CloseHandle( hProcess );
}

#elif __linux__

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

double CPUTimer::get_clock() const { 
  struct timespec et;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &et);
  return double(et.tv_sec) + 1.0e-9*double(et.tv_nsec);
}

double memory_usage(int peak){
  FILE* file = fopen("/proc/self/status", "r");
  char line[300], query[30], *match;

  if(peak == 1) strcpy(query, "VmPeak:");
  else strcpy(query, "VmSize:");
  
  if(file == NULL)
	  return -1.0;

  while (fgets(line, 300, file) != NULL)
    if ((match = strstr(line, query)) != NULL)
      return 1e3 * strtod(match+strlen(query), NULL);

  fclose(file);
  return -1.0;
}

#endif
