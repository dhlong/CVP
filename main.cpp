#include "cvp.h"

fstream iteration_report; // globally-accessed iteration report
fstream solve_report;     // globally-accessed solve report
SettingMapper settings;   // globally-accessed setting mapper

int main(int argc, char* argv[]){
  string ininame,inputname, outputname, timername, pajekname;
  Timer *timer = new CPUTimer;
  
  // get the system date and time
  // append date and time to the end of the report file name
  // so that the reports won't be overwritten after each run
  time_t rawtime;
  struct tm * timeinfo;
  char timestr[100];
  time(&rawtime );
  timeinfo = localtime (&rawtime);
  sprintf(timestr,"%04d%02d%02d%02d%02d",
	  timeinfo->tm_year+1900, 
	  timeinfo->tm_mon+1, 
	  timeinfo->tm_mday, 
	  timeinfo->tm_hour, 
	  timeinfo->tm_min);

  // start the timer
  timer->record();

  // Reading settings from file "CVP.ini"
  ininame = "CVP.ini";
  fstream fini(ininame.c_str(), fstream::in);
  settings.read(fini);
  fini.close();
  
  // generate input and report file names
  inputname = settings.gets("input file");
  outputname = inputname + "_iterations_" + string(timestr) + ".txt";
  timername  = inputname + "_solves_" + string(timestr) + ".txt";
  pajekname  = inputname + "pajek.net";
 
  // prepare files for reporting
  iteration_report.open(outputname.c_str(), fstream::out);
  solve_report.open(timername.c_str(), fstream::out);
  
  iteration_report<<"Solving problem \""<<inputname<<"\""<<endl<<endl;

  // Report settings
  iteration_report<<"Settings:"<<endl;
  settings.report(iteration_report);
  
  // check format of the input based on the filename given
  // if file name contains "grid" or "planar" --> genflot format
  // else tntp format
  FileFormat format = GENFLOT;
  if( (strstr(inputname.c_str(),"grid")   == NULL ) &&
      (strstr(inputname.c_str(),"planar") == NULL) )
    format = TNTP;
  
  MultiCommoNetwork net(inputname.c_str(), format);
  timer->record();
  
  cout << "Modelling" << endl;
  CVP_MCNF *aCVP =  new CVP_MCNF(new BPRFunction(net), net);
  //CVP_MCNF_KL *aCVP =  new CVP_MCNF_KL(net);

  timer->record();
  
  cout<<"Sovling"<<endl;
  aCVP->optimize();
  timer->record();

  cout<<"Finishing"<<endl;
  delete aCVP;
  timer->record();

  solve_report<<"Reading time    = "<<timer->elapsed(-4,-5)<<"s"<<endl;
  solve_report<<"Modelling time  = "<<timer->elapsed(-3,-4)<<"s"<<endl;
  solve_report<<"Solving time    = "<<timer->elapsed(-2,-3)<<"s"<<endl;
  solve_report<<"Finalizing time = "<<timer->elapsed(-1,-2)<<"s"<<endl;

  iteration_report << "Solving time = "
		   << left << setw(10) << setprecision(5) << fixed
		   << timer->elapsed(-2,-3) << "s"
		   << endl;

  iteration_report << "Peak memory  = " 
		   << left << setw(10) << setprecision(2) << fixed
		   << memory_usage()/1e3 << "kB" <<endl;

  iteration_report.close();
  solve_report.close();
  delete timer;
  return 0;
}

