#include "network.h"
#include <Eigen/Dense>


int main(){
  Timer timer;
  timer.record();
  string ininame = "CVP.ini";
  fstream fini(ininame.c_str(), fstream::in);
  SettingMapper settings;
  settings.read(fini);
  fini.close();

  MultiCommoNetwork net(settings.gets("input file").c_str(), GENFLOT);

  int V = net.getNVertex(), A = net.arcs.size();
  MatrixXd m = MatrixXd::Zero(V-1,A);
  FOR(a, A){
    if(net.arcs[a].head != V-1) m(net.arcs[a].head, a) = 1.0;
    if(net.arcs[a].tail != V-1) m(net.arcs[a].tail, a) = -1.0;
  }
  
  timer.record();
  MatrixXd mmT = m*m.transpose();

  timer.record();
  MatrixXd inv = mmT.inverse();

  timer.record();
  MatrixXd I = mmT*inv;
  
  timer.record();
  MatrixXd M = m.transpose()*I*m;
  timer.record();

  double eps = 0.0;
  double averageM = 0.0;
  int count = 0;
  FOR(i, A) FOR(j, A) averageM += fabs(M(i,j));
  averageM /= (double(A)*A);
  FOR(i, A) FOR(j, A) M(i, j) /= averageM;
  FOR(i, A) FOR(j, A) if(fabs(M(i,j)) > 1e-6) count++;
  FOR(i, V-1) FOR(j, V-1) 
    if(i == j) updatemax(eps, fabs(I(i,j)-1));
    else updatemax(eps, fabs(I(i,j)));

  cout<<"eps max = "<<eps<<endl;
  cout<<"average M = "<<averageM<<endl;
  cout<<"Non zeros: "<<count<<"/"<<A*A<<endl;
  cout<<"Multiplication 2 time "<<setprecision(6)<<timer.elapsed()<<"s"<<endl;
  cout<<"Multiplication 1 time "<<setprecision(6)<<timer.elapsed(-4,-5)<<"s"<<endl;
  cout<<"Inversion time "<<setprecision(6)<<timer.elapsed(-4,-3)<<"s"<<endl;
}
