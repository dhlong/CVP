#include "cvp.h"

fstream iteration_report;
fstream solve_report;
SettingMapper settings;

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
	  timeinfo->tm_mon, 
	  timeinfo->tm_mday, 
	  timeinfo->tm_hour, 
	  timeinfo->tm_min);

  // settings information
  // if file CVP.ini does not exist, default settings will be used
  ininame = "CVP.ini";
  fstream fini(ininame.c_str(), fstream::in);
  settings.read(fini);
  fini.close();

  // start the timer
  timer->record();

  cout<<"Mark1"<<endl;

  inputname = settings.gets("input file");
  cout<<"Input file to be read: "<<inputname<<endl;
  outputname = inputname + "_iterations_" + string(timestr) + ".txt";
  timername  = inputname + "_solves_" + string(timestr) + ".txt";
  pajekname  = inputname + "pajek.net";
 
  iteration_report.open(outputname.c_str(), fstream::out);
  solve_report.open(timername.c_str(), fstream::out);
  
  cout<<"Solving problem \""<<inputname<<"\""<<endl<<endl;
  iteration_report<<"Solving problem \""<<inputname<<"\""<<endl<<endl;
    
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

  fini.open("CVP.ini", fstream::in);
  aCVP->read_settings(fini);
  fini.close();
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
		   << right << setw(10) << setprecision(5) << fixed
		   << timer->elapsed(-2,-3) << "s"
		   << endl;

  memory_report(iteration_report);

  iteration_report.close();
  solve_report.close();
  delete timer;
  return 0;
}
