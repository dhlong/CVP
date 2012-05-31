#include "cvp.h"

fstream iteration_report;
fstream solve_report;

int main(int argc, char* argv[]){
  string ininame,inputname, outputname, timername, pajekname;
  Timer timer;
  Real tic1, tic2, tic3, tic4, tic5, tic6;
  
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

  // start the timer
  tic1 = timer.elapsed();

  read_param(ininame, "input file", inputname);
  outputname = inputname + "_iterations_" + string(timestr) + ".txt";
  timername  = inputname + "_solves_" + string(timestr) + ".txt";
  pajekname  = inputname + "pajek.net";
 
  iteration_report.open(outputname.c_str(), fstream::out);
  solve_report.open(timername.c_str(), fstream::out);
  

  iteration_report<<"Solving problem "<<inputname<<endl<<endl;
    
  // check format of the input based on the filename given
  // if file name contains "grid" or "planar" --> genflot format
  // else tntp format
  FileFormat format = GENFLOT;
  if( (strstr(inputname.c_str(),"grid")   == NULL ) &&
      (strstr(inputname.c_str(),"planar") == NULL) )
    format = TNTP;
  
  MultiCommoNetwork net(inputname.c_str(), format);
  tic2 = timer.elapsed();
  
  cout<<"Modelling"<<endl;
  //CVP_MCNF *aCVP =  new CVP_MCNF(new BPRFunction(net), net);
  CVP_MCNF_KL *aCVP =  new CVP_MCNF_KL(net);

  fstream fini("CVP.ini", fstream::in);
  aCVP->read_settings(fini);
  fini.close();
  tic3 = timer.elapsed();
  
  cout<<"Sovling"<<endl;
  aCVP->optimize();
  tic4 = timer.elapsed();

  cout<<"Finishing"<<endl;
  delete aCVP;

  tic5 = timer.elapsed();

  solve_report<<"Reading time: "<<tic2-tic1<<endl;
  solve_report<<"Modelling time:\t" <<tic3-tic2<<endl;
  solve_report<<"Solving time:\t"   <<tic4-tic3<<endl;
  solve_report<<"Finalizing time:\t"<<tic5-tic4<<endl;

  iteration_report << "Solving time     : "
		   << right << setw(10) << setprecision(5) << fixed << tic4-tic3<<"s"
		   << endl;

  memory_report(iteration_report);

  iteration_report.close();
  solve_report.close();
  return 0;
}
