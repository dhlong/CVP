
#include "cvp.h"

fstream iteration_report;
fstream solve_report;

int main(int argc, char* argv[]){
  string outputname, timername, pajekname;
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
  fstream fini("CVP.ini", fstream::in);

  // start the timer
  tic1 = timer.elapsed();




  /**
   * set the reports' file names
   * if no argument is given, the 400x400 grid network will be generated
   * if more than 1 arguments are given:
   * +  the 1st argument will be the path to the input file
   *        (in case of tntp format, it is the path excluding 
   *        the "_net.txt" and "_trips.txt" parts of the file names)
   * +  the 2nd argument will be the file name of the iteration report
   * +  the 3rd argument will be the file name of the iteration report   
   * If the 2nd and 3rd arguments are not specified, they will be assigned
   * based on the 1st argument
   *
   */

  if(argc == 1)
    outputname = "CVPiterations.txt",
      timername = "CVPsolves.txt",
      pajekname = "CVPpajek.net";
  else {
    if(argc>=3) outputname = string(argv[2]);
    else outputname = string(argv[1]) + "_iterations_" + string(timestr) + ".txt";

    if(argc>=4) timername = string(argv[3]);
    else timername = string(argv[1]) + "_solves_" + string(timestr) + ".txt";

    pajekname = string(argv[1]) + "pajek.net";
  }
 
  iteration_report.open(outputname.c_str(), fstream::out);
  solve_report.open(timername.c_str(), fstream::out);
  

  // If no argument is specified
  // solve the 400x400 generated grid network
  if(argc == 1) {
    Network net(400, 10.0);
    net.write_pajek(pajekname.c_str());
    tic2 = timer.elapsed();

    CVP *aCVP =  new CVP_NF(new QuarticFunction(net), net);
    aCVP->read_settings(fini);
    tic3 = timer.elapsed();

    aCVP->optimize();
    tic4 = timer.elapsed();

    delete aCVP;
    tic5 = timer.elapsed();
  } 
  else {
    iteration_report<<"Solving problem "<<argv[1]<<endl;
    
    // check format of the input based on the filename given
    // if file name contains "grid" or "planar" --> genflot format
    // else tntp format
    FileFormat format = GENFLOT;
    if(strstr(argv[1],"grid") == NULL && strstr(argv[1],"planar") == NULL)
      format = TNTP;

    MultiCommoNetwork net(argv[1], format);

    cout<<"Generating pajek"<<endl;
    net.write_pajek(pajekname.c_str());
    tic2 = timer.elapsed();

    cout<<"Modelling"<<endl;
    CVP_MCNF *aCVP =  new CVP_MCNF(new BPRFunction(net), net);
    //CVP_MCNF *aCVP =  new CVP_MCNF(new KleinrockFunction(net), net);
    aCVP->read_settings(fini);
    tic3 = timer.elapsed();

    cout<<"Sovling"<<endl;
    aCVP->optimize();
    //aCVP->solve_by_dijkstra();
    tic4 = timer.elapsed();

    cout<<"Finishing"<<endl;
    delete aCVP;
    tic5 = timer.elapsed();
  }

  solve_report<<"Reading time: "<<tic2-tic1<<endl;
  solve_report<<"Modelling time:\t" <<tic3-tic2<<endl;
  solve_report<<"Solving time:\t"   <<tic4-tic3<<endl;
  solve_report<<"Finalizing time:\t"<<tic5-tic4<<endl;

  /*
  iteration_report << "Reading time     : "
		   << right << setw(10) << setprecision(5) << fixed << tic2-tic1
		   << "s"<<endl;
  iteration_report << "Modelling time   : "
		   << right << setw(10) << setprecision(5) << fixed << tic3-tic2
		   <<"s"<<endl;
  */
  iteration_report << "Solving time     : "
		   << right << setw(10) << setprecision(5) << fixed << tic4-tic3<<"s"
		   <<endl;

  memory_report(iteration_report);

  iteration_report.close();
  solve_report.close();
  return 0;
}
