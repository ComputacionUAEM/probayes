/*=============================================================================
 * Product        : ProBT 
 * File           : RedBayesiana.cpp
 * Author         : Esau Escobar
 * Creation       : Lun Abr 23 14:41:15 2012
 *
 *=============================================================================
 *        (c) Copyright 2012  -  all rights reserved
 *=============================================================================
*/

#include <pl.h>
#include <iostream>
#include <math.h>
using namespace std;

//Funciones externas
float recta(float x,float m,float dm,float M,float dM){
  float y,p;
  p= (M-m)/(dM-dm);
  y=(p*x)+(m-(p*dm));
  return y;
}

void f_mean(plValues &mean, const plValues &n){
  if(n[0]-n[1]<5&&n[0]-n[1]>-5)
    mean[0]=0;
  else if(n[0]-n[1]>20)
    mean[0]=9;
  else if(n[0]-n[1]<-20)
    mean[0]=-9;
  else if(n[0]-n[1]>0)
    mean[0]=recta(n[0]-n[1],1,5,9,20);
  else
    mean[0]=recta(n[0]-n[1],-1,-5,-9,-20);
}

void f_std(plValues &std, const plValues &n){
  float d = abs((float)(n[0]-n[1]));
  cout<<"d:"<<d<<endl;
  if(d<5)
    std[0]=0.1;
  else if(d>20)
    std[0]=0.8;
  else
    std[0]=recta(d,0.1,5,2.0,20);
  cout<<"std:"<<std[0]<<endl;
}

int main() 
{
  /**********************************************************************
   VARIABLES SPECIFICATION
  **********************************************************************/
  plIntegerType t_rank(40,150);      //El rango de Z es de 40cm a 1.5mts
  plIntegerType vt_rank(-9,9);       //28 posibles combinaciones de
                                     //  Velocidad * Tiempo
  plSymbol sent("St",t_rank); //Espacio de Situacion Sensorial t
  plSymbol sent1("St+1",t_rank); //Espacio de Situacion Sensorial t
  plSymbol vt("VT",vt_rank);

  plValues result(sent^sent1^vt);
  plVariable sens(sent^sent1);
  /**********************************************************************
   PARAMETRIC FORM SPECIFICATION
  **********************************************************************/
  //P(sent) = uniform(sent)
  plUniform p_sent(sent);
  plUniform p_sent1(sent1);

  //Funciones externas para computar la media y la desv.est.
  plExternalFunction f_mu(sens,f_mean);
  plExternalFunction f_sd(sens,f_std);

  //plExternalFunction f_mu1(sent,f_mean1);
  //plExternalFunction f_sd1(sent,f_std1);

  //P(vt | t,t1)=CndBellShape(t,t1,f_mu(st,st1),f_sd(st,st1))
  //plCndBellShape p_sent1(sent1,sent,f_mu1,f_sd1);
  plCndBellShape p_vt(vt,sent^sent1,f_mu,f_sd);

  // cerr<<"t_rank "<<<<"\n";
  // cerr<<"vt_rank "<<<<"\n";
  // cerr<<"sent "<<<<"\n";
  // cerr<<"sent1 "<<<<"\n";
  // cerr<<"vt "<<<<"\n";
  
  /**********************************************************************
   DECOMPOSITION
  **********************************************************************/
  //P(vt sent sent1) = P(vt | sent sent1)P(sent)P(sent1)
  plJointDistribution jd(sent^sent1^vt, p_sent*p_sent1*p_vt);
  jd.draw_graph("red_bayesiana.fig");
  cout<<"OK\n";
  /**********************************************************************
   PROGRAM QUESTION
  **********************************************************************/
  plCndDistribution cnd_question;
  plDistribution question;
  int st,v,st1;

  
  //Dame P(sent1| vt sent)
  jd.ask(cnd_question,sent1,sent^vt);
  cout<<"dame sent: "; cin>>st;
  result[sent] = st;
  cout<<"\ndame st1: "; cin>>st1;
  result[sent1] = st1;
  
  
  plDistribution compiled;
  p_vt.instantiate(question,result);
  question.compile(compiled);
  cout<<"Dist p_vt compilada:\n"<<compiled<<endl;
  compiled.plot("archivo",PL_DEFAULT_PLOT,19);
  
  cout<<"\ndame vt: "; cin>>v;
  result.reset();
  result[sent] = st;
  result[vt] = v;
  //Dame P(sent1|vt=v sent=st)
  cnd_question.instantiate(question,result);
  plDistribution compiled_question;
  cout<<"La dist instanciada:\n"<<question<<endl;


  question.compile(compiled_question);
  cout<<"La distribucion es :\n"<<compiled_question<<endl;
  compiled_question.best(result);
  compiled_question.plot("archivo2",PL_DEFAULT_PLOT,110);
  cout<<"\nLo mas probable es: "<<result[sent1]<<"\n";
  
  return 0;
}
