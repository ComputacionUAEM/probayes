/*=============================================================================
 * Product        : ProBT 
 * File           : RedBayesiana.cpp
 * Author         : Esau Escobar
 * Creation       : Lun Abr 25 14:00:00 2012
 *
 *=============================================================================
 *        (c) Copyright 2012  -  all rights reserved
 *=============================================================================
*/

#include <pl.h>
#include <iostream>
#include <math.h>

using namespace std;

void f_mean1(plValues &mean, const plValues &n){
  float ang = n[0];
  mean[0]= (n[1]*sin(ang))+n[2];
}

void f_std1(plValues &std, const plValues &n){
  std[0]=2.0;
}

void f_mean2(plValues &mean, const plValues &n){
  float ang = n[0];
  mean[0]= (n[1]*cos(ang))+n[2];
}

void f_std2(plValues &std, const plValues &n){
  std[0]=2.0;
}

void f_mean3(plValues &muthe, const plValues &n){
   muthe[0]= n[1]+n[0];
}

void f_std3(plValues &std, const plValues &n){
  std[0]=2.0;
}

int main() 
{
  /**********************************************************************
   VARIABLES SPECIFICATION
  **********************************************************************/
  plIntegerType x_rank(-77,77);     //El rango de X es de -77cm a 77cm
  plIntegerType z_rank(0,200);      //El rango de Z es de 0cm a 2mts
  plRealType vt_rank(-18.0,18.0,500);     //Rango de velocidad traslacional
  plRealType vr_rank(-1.26,1.26,500);     //Rango vel. rot. en rad
  plRealType theta_rank(-6.28,6.28,500);  //Ranog de theta

  plSymbol sx("x",x_rank);     //Variable de X tiempo t
  plSymbol sx1("x1",x_rank);   //Variable de X tiempo t+1
  plSymbol sz("z",z_rank);     //Variable de Z tiempo t
  plSymbol sz1("z1",z_rank);   //Variable de Z tiempo t+1
  plSymbol vr("Vrot",vr_rank);
  plSymbol vt("Vtras",vt_rank);
  plSymbol theta("Theta",theta_rank);
  plSymbol theta1("Theta1",theta_rank);

  //plValues result(sent^sent1^vt);
  plVariable stx(theta1^vt^sx);
  plVariable stz(theta1^vt^sz);
  plVariable sth(theta^vr);

  plValues result(sx^sx1^sz^sz1^vt);
  plValues result2(theta^theta1^vr);
  plValues resulth(theta^vr^theta1);

  plValues results(sx^sx1^theta1^vt);
  plValues resultb(sz^sz1^theta1^vt);
  /**********************************************************************
   PARAMETRIC FORM SPECIFICATION
  **********************************************************************/
  //P(sx) = uniform(sx)
  plUniform p_sx(sx);
  //P(sz) = uniform(sz)
  plUniform p_sz(sz);
  //P(vt) = uniform(vt)
  plCUniform p_vt(vt);
  //P(vr) = uniform(vr)
  plCUniform p_vr(vr);
  //P(theta)=uniform(theta)
  plCUniform p_theta(theta);

  plExternalFunction f_muth(sth,f_mean3);
  plExternalFunction f_sdth(sth,f_std3);

  //P(theta1 | vr,theta)=CndBellShape(theta,theta^vr,f_mu(vr),f_sd(vr))
  plCndNormal p_theta1(theta1,theta^vr,f_muth,f_sdth);

  //Funciones externas para computar la media y la desv.est.
  plExternalFunction f_mux(stx,f_mean1);
  plExternalFunction f_sdx(stx,f_std1);

  //P(sx1 | theta,vt,x)=CndBellShape(sx1,th^vt^x,f_mu(th,vt,sx),f_sd(th,vt,sx))
  plCndBellShape p_sx1(sx1,theta1^vt^sx,f_mux,f_sdx);

  //Funciones externas para computar la media y la desv.est.
  plExternalFunction f_muz(stz,f_mean2);
  plExternalFunction f_sdz(stz,f_std2);

  //P(sz1 | theta,vt,z)=CndBellShape(sz1,th^vt^z,f_mu(th,vt,sz),f_sd(th,vt,sz))
  plCndBellShape p_sz1(sz1,theta1^vt^sz,f_muz,f_sdz);

  /**********************************************************************
   DECOMPOSITION
  **********************************************************************/
  //P(sx sz sx1 sz1 vt vr) = P(sx)P(sz)P(vt)P(vr)P(sx1|vr,vt,sx)P(sz1|vr,vt,sz)
  plJointDistribution jd(sx^sz^vr^vt^theta^theta1^sx1^sz1, p_sx*p_sz*p_vr*p_vt*p_theta*p_theta1*p_sx1*p_sz1);
  jd.draw_graph("red_bayesiana.fig");
  cout<<"OK\n";
  /**********************************************************************
   PROGRAM QUESTION
  **********************************************************************/
  int x,x1,z,z1,mc=200;
  float vr2,vt2,theta2;
  
  cout<<"Prueba de MD";
  cout<<"\ndame sx: "; cin>>x;
  results[sx] = x; 
  cout<<"\ndame sz: "; cin>>z;
  resultb[sz] = z;
  cout<<"\ndame theta: "; cin>>theta2;
  resulth[theta] = theta2;
  cout<<"\ndame vr: "; cin>>vr2;
  resulth[vr] = vr2;
  cout<<"\ndame vt: "; cin>>vt2;
  results[vt] = vt2;
  resultb[vt] = vt2;
  
  plDistribution compiled6,question6;
  p_theta1.instantiate(question6,resulth);
  question6.compile(compiled6);
  //cout<<"Dist p_theta1 compilada:\n"<<compiled6<<endl;
  compiled6.best(results);
  compiled6.best(resultb);
  
  plDistribution compiled,question1;
  p_sx1.instantiate(question1,results);
  question1.compile(compiled);
  //cout<<"Dist p_sx1 compilada:\n"<<compiled<<endl;
  compiled.best(results);
  cout<<"\nSx1 mas probable es: "<<results[sx1]<<"\n";
  
  plDistribution compiled2,question2;
  p_sz1.instantiate(question2,resultb);
  question2.compile(compiled2);
  //cout<<"Dist p_sz1 compilada:\n"<<compiled2<<endl;
  compiled2.best(resultb);
  cout<<"Sz1 mas probable es: "<<resultb[sz1]<<"\n";

  cout<<"\nArmar MI";
  cout<<"\ndame sx: "; cin>>x;
  result[sx] = x;
  //result2[sx] = x;
  cout<<"\ndame sx1: "; cin>>x1;
  result[sx1] = x1;
  //result2[sx1] = x1;
  cout<<"\ndame sz: "; cin>>z;
  result[sz] = z;
  //result2[sz] = z;
  cout<<"\ndame sz1: "; cin>>z1;
  result[sz1] = z1;
  //result2[sz1] = z1;
  cout<<"\ndame theta: "; cin>>theta2;
  result2[theta] = theta2;
  cout<<"\ndame theta1: "; cin>>theta2;
  result2[theta1] = theta2;
  //result2[theta] = theta2;
  //Dame P(vt| x z x1 z1)
  plDistribution compiled3,question3;
  plCndDistribution cnd_question;
  jd.ask(cnd_question,vt,sx^sz^sx1^sz1,mc);
  cnd_question.instantiate(question3,result);
  question3.n_compile(compiled3,200);
  //cout<<"Dist p_vt compilada:\n"<<compiled3<<endl;
  compiled3.plot("Vtras",PL_DEFAULT_PLOT,200);
  compiled3.best(result);
  cout<<"\nVt mas probable es: "<<result[vt]<<"\n";
  
  //Dame P(vr| theta, theta1)
  plCndDistribution cnd_question2;
  plDistribution compiled4,question4;
  jd.ask(cnd_question2,vr,theta^theta1,mc);
  cnd_question2.instantiate(question4,result2);
  question4.n_compile(compiled4,200);
  //cout<<"Dist p_vt compilada:\n"<<compiled4<<endl;
  compiled4.plot("Vrot",PL_DEFAULT_PLOT,200);
  compiled4.best(result2);
  cout<<"\nVr mas probable es: "<<result2[vr]<<"\n";
  
  return 0;
}
