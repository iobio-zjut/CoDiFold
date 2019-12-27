//@pengchunxiang


#include "mymethod.hh"

#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <cmath>
#include <string.h>


#include <protocols/abinitio/ClassicAbinitio.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/rms_util.hh>


void Mymethod::computed_coefficient_contact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double **contact1,double **contact2){
    
     double temp1=0.0;
     double temp2=0.0;
     double temp3=0.0;
     for(int i = 0 ; i < int(pd_position_first.size()) ; i++){
	temp1=(contact1[pd_position_first[i]][pd_position_second[i]]+contact2[pd_position_first[i]][pd_position_second[i]])/2.0;
// 	std::cout<<contact1[pd_position_first[i]][pd_position_second[i]]<<contact2[pd_position_first[i]][pd_position_second[i]]<<std::endl;
	temp2=fabs(contact1[pd_position_first[i]][pd_position_second[i]]-contact2[pd_position_first[i]][pd_position_second[i]]);
	temp3=temp1/exp(temp2);
	coefficient_contact.push_back(temp3);
	
    }
  
}

void Mymethod::energy_distanceprofile(std::vector<float>& coefficient_contact,std::vector<float>& decoy_ditance,std::vector <float>& pd_peak_distance){
      float energy_dp=0.0,sum=0.0;
      sum_energy_distanceprofile=0.0;
      std::ofstream outfile("/home/pcx/dp.txt");
  /*    if(!outfile.is_open()){
      
	std::cout<<"open failed";
	outfile.close();
      }*/
      for(int i = 0 ;i< int(coefficient_contact.size()) ; i++){
      
	energy_dp=fabs(decoy_ditance[i]-pd_peak_distance[i]);
	std::cout<<"energy_dp=   "<<energy_dp<<"  coefficient_contact=   "<<coefficient_contact[i]<<std::endl;
	sum_energy_distanceprofile=sum_energy_distanceprofile+coefficient_contact[i]*energy_dp;
	sum=sum+energy_dp;
	outfile<<decoy_ditance[i]<<" "<<pd_peak_distance[i]<<" "<<std::endl;
      }
      outfile.clear();
      outfile.close();
      std::cout<<"average_err= "<<sum/(int(decoy_ditance.size()))<<"sum=  "<<sum<<std::endl;
      
}

void Mymethod::computed_bothexist_dpandcontact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double contact1[][108],double contact2[][108]){
  
    for(int i = 0 ; i < int(pd_position_first.size()) ; i++){
    
      if(contact1[pd_position_first[i]][pd_position_second[i]]!=0.0&&contact2[pd_position_first[i]][pd_position_second[i]]!=0.0)
	both_exist_number++;
    }
  }
  
  
  
  void Mymethod::computed_distanceprofile(std::vector<float>& decoy_ditance,std::vector <float>& pd_peak_distance){
      float energy_dp=0.0,sum=0.0;
      std::ofstream outfile("/home/pcx/protein/1ALY/dp.txt");
  /*    if(!outfile.is_open()){
      
	std::cout<<"open failed";
	outfile.close();
      }*/
      for(int i = 0 ;i< int(decoy_ditance.size()) ; i++){
      
	energy_dp=fabs(decoy_ditance[i]-pd_peak_distance[i]);
	
	sum=sum+energy_dp;
	
	outfile<<decoy_ditance[i]<<" "<<pd_peak_distance[i]<<" "<<std::endl;
      }
      outfile.clear();
      outfile.close();
      std::cout<<"average_err= "<<sum/(int(decoy_ditance.size()))<<"sum=  "<<sum<<std::endl;
      
}

 void Mymethod::show_dppair_contact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double contact1[][108],double contact2[][108])
{
  
  std::ofstream outfile("/home/pcx/contact/2LZMA/dppair_contact.txt");
  for(int i = 0 ;i<int(pd_position_first.size()) ; i++){
  
    outfile<<contact1[pd_position_first[i]+1][pd_position_second[i]+1]<<" "<<contact2[pd_position_first[i]+1][pd_position_second[i]+1]<<" "<<std::endl;
  }
  outfile.clear();
  outfile.close();
}


void Mymethod::dp_initialization(core::pose::Pose& pose,core::pose::Pose mypose[] )
{	
  bool success=false;
  
  for(int i = 0 ;i<pop_size ; i++){
    
 //   success = do_stage1_cycles( pose);
    if(success){
      std::cout<<i<<"     fragment assemble success"<<std::endl;
    mypose[i]=pose;
    }else
      std::cout<<i<<"     fragment assemble failed"<<std::endl;
  }

}
void Mymethod::dp_mutation(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand3, rand_fra1, rand_fra2;
  float rand_cr;
  core::pose::Pose pose1,pose2;
  bool success = false;
  bool exchange = false;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  
  rand_fra1 = rand()%(seq_length-9)+1;
  do rand_fra2 = rand()%(seq_length-9)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  posemutation = mypose[rand3];
  for(int k = 0 ; k<9 ; k++){
    posemutation.set_phi(rand_fra1 + k,pose1.phi(rand_fra1 + k));
    posemutation.set_psi(rand_fra1 + k,pose1.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<9 ; k++){
    posemutation.set_phi(rand_fra2 + k,pose2.phi(rand_fra2 + k));
    posemutation.set_psi(rand_fra2 + k,pose2.psi(rand_fra2 + k));
  }
}


void Mymethod::dp_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length){
  int cr_pos;
  cr_pos = rand()%(seq_length-9)+1;
    
  for(int k = 0 ; k<9 ; k++){
    posemutation.set_phi(cr_pos + k, target.phi(cr_pos + k));
    posemutation.set_psi(cr_pos + k, target.psi(cr_pos + k));
  }
}
 
    
bool Mymethod::average_distance_error_discriminate(const core::pose::Pose& pose_trial,const core::pose::Pose& pose_target,std::vector<int>& pd_position_first,
			    std::vector<int>& pd_position_second,std::vector<float>& pd_peak_distance,std::vector<float>& decoy_ditance_trial,std::vector<float>& decoy_ditance_target){
  
  std::vector<float> delta_trial,delta_target;
  float final_delta_trial=0.0,final_delta_target=0.0;
 

  
  for(int k = 0;k < int(pd_peak_distance.size());k++){
    delta_trial.push_back(fabs(decoy_ditance_trial[k]-pd_peak_distance[k]));
    if(delta_trial[k]>6.5)
      delta_trial[k]=6.5;
//      std::cout<<"delta_trial    "<<delta_trial[k]<<std::endl;
  }
  
  for(int k = 0;k < int(pd_peak_distance.size());k++){
    delta_target.push_back(fabs(decoy_ditance_target[k]-pd_peak_distance[k]));
    if(delta_target[k]>6.5)
      delta_target[k]=6.5;
//      std::cout<<"delta_target    "<<delta_target[k]<<std::endl;    
  }

  for(int k = 0 ;k<int(delta_trial.size());k++)
  final_delta_trial=final_delta_trial+delta_trial[k];
  
  final_delta_trial=final_delta_trial/(int(delta_trial.size()));
  
  for(int k = 0 ;k<int(delta_target.size());k++ )
  final_delta_target=final_delta_target+delta_target[k];
  
  final_delta_target=final_delta_target/(int(delta_target.size()));
 // std::cout<<"final_delta_trial    "<<final_delta_trial<<std::endl;
 // std::cout<<"final_delta_target    "<<final_delta_target<<std::endl;
  if(final_delta_trial<final_delta_target){
    std::vector<float> pk;
    float paccept=0.0,rand5;
    
    
    for(int k = 0 ; k < int(delta_trial.size()) ; k++){
    
     // pk.push_back(exp(-delta_trial[k]/belta));
      pk.push_back(exp(-delta_trial[k]/4));
      std::cout<<"pk    "<<pk[k]<<std::endl;
    }
    for(int k = 0 ;k < int(pk.size());k++){
    
      paccept=pk[k]+paccept;
  //    std::cout<<"paccept=pk[k]+paccept    "<<paccept<<"pk    "<<pk[k]<<std::endl;
    }
    paccept=paccept/int(pk.size());
    std::vector<float>().swap(pk);  
    std::vector<float>().swap(delta_trial);
    std::vector<float>().swap(delta_target);    
 //   std::cout<<"paccept!!!!!"<<paccept<<std::endl;
    rand5=rand()%100/(double)101;
    if(rand5<paccept){
    std::cout<<"exchange    !!!!"<<std::endl;
      return true;
      
    }else{
      std::cout<<"retain target!!!"<<std::endl;
      return false;
    }
    
  }else{
    
    std::cout<<"average distance error retain target!!!"<<std::endl;
  
    return false;
  }
  
  
  
  
}
void Mymethod::energy_distanceprofile2(std::vector<float>& decoy_ditance,std::vector<float>& coefficient_contact,float& sumE,std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925){
  float temp=0.0;
  sumE=0.0;
      for(int i = 0 ;i<int(decoy_ditance.size()); i++){
      if(decoy_ditance[i]<8){
	temp=-distance_number375[i]/fabs(decoy_ditance[i]-3.75)-distance_number425[i]/fabs(decoy_ditance[i]-4.25)-distance_number475[i]/fabs(decoy_ditance[i]-4.75)-distance_number525[i]/fabs(decoy_ditance[i]-5.25)
	-distance_number575[i]/fabs(decoy_ditance[i]-5.75)-distance_number625[i]/fabs(decoy_ditance[i]-6.25)-distance_number675[i]/fabs(decoy_ditance[i]-6.75)
	-distance_number725[i]/fabs(decoy_ditance[i]-7.25)-distance_number775[i]/fabs(decoy_ditance[i]-7.75);
      }else{
      
	temp=log(decoy_ditance[i]-8);
      }
	sumE=sumE+coefficient_contact[i]*temp;
      }
  std::cout<<"sumE     "<<sumE<<std::endl;
}

/////////////////////////////////////////////////////////////////////////

void Mymethod::energy_distanceprofile3(std::vector<float>& decoy_ditance,std::vector<double>& coefficient_contact,double& sumE, std::vector<int>& pd_totalNumber,std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925){
  float temp=0.0;
  sumE=0.0;
      for(int i = 0 ;i<int(decoy_ditance.size()); i++){
	temp=-distance_number375[i]/exp(pow((decoy_ditance[i]-3.75),2)) - distance_number425[i]/exp(pow((decoy_ditance[i]-4.25),2)) - distance_number475[i]/exp(pow((decoy_ditance[i]-4.75),2))-distance_number525[i]/exp(pow((decoy_ditance[i]-5.25),2))
	-distance_number575[i]/exp(pow((decoy_ditance[i]-5.75),2)) - distance_number625[i]/exp(pow((decoy_ditance[i]-6.25),2)) - distance_number675[i]/exp(pow((decoy_ditance[i]-6.75),2))
	-distance_number725[i]/exp(pow((decoy_ditance[i]-7.25),2)) - distance_number775[i]/exp(pow((decoy_ditance[i]-7.75),2)) - distance_number825[i]/exp(pow((decoy_ditance[i]-8.25),2)) - distance_number875[i]/exp(pow((decoy_ditance[i]-8.75),2))
	-distance_number925[i]/exp(pow((decoy_ditance[i]-9.25),2));
     
	//std::cout<<"pd_totalNumber["<<i<<"]="<<pd_totalNumber[i]<<"   decoy_ditance size =      "<<decoy_ditance.size()<<"    pd_totalNumber size =      "<<pd_totalNumber.size()<<std::endl;
	sumE += (coefficient_contact[i]/pd_totalNumber[i])*temp;
	//sumE += ( coefficient_contact[i]) * temp;
      }
     // sumE = sumE / double(decoy_ditance.size());
 // std::cout<<"sumE     "<<sumE<<std::endl;
}



///////////////////////////////////////////////////////////////////////

void Mymethod::energy3(std::vector<float>& decoy_ditance, std::vector<int>& pd_totalNumber, float& sumE, std::vector<float>& E_res, std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925){
    float temp = 0.0;
    sumE = 0.0;
    for(int i = 0; i < int(decoy_ditance.size()); i++){
	temp = -distance_number375[i]/fabs(decoy_ditance[i]-3.75) - distance_number425[i]/fabs(decoy_ditance[i]-4.25) - distance_number475[i]/fabs(decoy_ditance[i]-4.75) - distance_number525[i]/fabs(decoy_ditance[i]-5.25)
	-distance_number575[i]/fabs(decoy_ditance[i]-5.75) - distance_number625[i]/fabs(decoy_ditance[i]-6.25) - distance_number675[i]/fabs(decoy_ditance[i]-6.75)
	- distance_number725[i]/fabs(decoy_ditance[i]-7.25) - distance_number775[i]/fabs(decoy_ditance[i]-7.75) - distance_number825[i]/fabs(decoy_ditance[i]-8.25)
	- distance_number875[i]/fabs(decoy_ditance[i]-8.75);

	sumE = sumE + temp/pd_totalNumber[i];
	E_res.push_back(temp/pd_totalNumber[i]);
    }
    std::cout<<"res_size:   "<<int(E_res.size()) <<std::endl;
    sumE = sumE/int(pd_totalNumber.size());
    std::cout<<"sumE:   "<<sumE<<int(E_res.size()) <<std::endl;
}

bool Mymethod::dp_contact_select(std::vector<float>& dp_E_trial, std::vector<float>& dp_E_target, float sumE_trial, float sumE_target){
    double p_accept = 0.0;
    int beta = 4;
    if(sumE_trial < sumE_target){
	double temp = 0.0;
	for(int k = 0; k < int(dp_E_trial.size()); k++){
	    temp = exp(-fabs(dp_E_trial[k] - dp_E_target[k])/beta);
	    p_accept = p_accept+temp;
	}
	p_accept = p_accept/int(dp_E_trial.size());

	if(rand()%100/(double)101 < p_accept)
	    return true;
	else 
	    return false;
    }
    else 
	return false;
}

void Mymethod::loop_crossover(core::pose::Pose mypose[], core::pose::Pose& posemutation, std::vector<int>& loop_pos_start, std::vector<int>& loop_pos_end, int pop_size, int target_id){
    int rand1, rand2;
    do rand1 = rand()%pop_size; while(rand1 == target_id);
    rand2 = rand()%int(loop_pos_start.size());
    posemutation = mypose[target_id];
    
    for(int i = loop_pos_start[rand2]; i <= loop_pos_end[rand2]; i++){
      posemutation.set_phi(i+1, mypose[rand1].phi(i+1));
      posemutation.set_psi(i+1, mypose[rand1].phi(i+1));
      posemutation.set_omega(i+1, mypose[rand1].phi(i+1));
    }
  
}



void Mymethod::set_default()
{
  int pop_size = 50;
  int belta=4;
  float CR=0.5;
}


///////////convex
double Mymethod::find_the_closest(double trial[], double decoy[][600], int seq_length, int kk)
{
  double temp= 0.0;
  for(int i = 0; i < seq_length * 3 ; i++){
  
    temp = temp + (trial[i]-decoy[kk][i]) * (trial[i]-decoy[kk][i]);
  }
  return std::sqrt(temp);
}

double Mymethod::find_the_closest_zxg(double* trial, double** decoy, int seq_length, int kk)
{
  double temp= 0.0;
  for(int i = 0; i < seq_length * 3 + 1; i++){
  
    temp = temp + (trial[i]-decoy[kk][i]) * (trial[i]-decoy[kk][i]);
  }
  return std::sqrt(temp);
}


void Mymethod::convex_mutation1(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand3, rand_fra1, rand_fra2;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  pose3 = mypose[rand3];
  for(int k = 0 ; k<3 ; k++){
    pose3.set_phi(rand_fra1 + k,pose1.phi(rand_fra1 + k));
    pose3.set_psi(rand_fra1 + k,pose1.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose3.set_phi(rand_fra2 + k,pose2.phi(rand_fra2 + k));
    pose3.set_psi(rand_fra2 + k,pose2.psi(rand_fra2 + k));
  }
  posemutation = pose3;
}

void Mymethod::convex_mutation2(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand3, rand4, rand_fra1, rand_fra2, rand_fra3;;
  core::pose::Pose pose;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  do rand4 = rand()%pop_size; while(rand4 ==rand3 || rand4 == rand2 || rand4 == rand1 || rand4 == target_id);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1;while(rand_fra3 ==rand_fra2 || rand_fra3 == rand_fra1);
   
  pose = mypose[rand4];
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    pose.set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra2 + k, mypose[rand2].phi(rand_fra2 + k));
    pose.set_psi(rand_fra2 + k, mypose[rand2].psi(rand_fra2 + k));
  }
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra3 + k, mypose[rand3].phi(rand_fra3 + k));
    pose.set_psi(rand_fra3 + k, mypose[rand3].psi(rand_fra3 + k));
  }  
  posemutation = pose;
}

void Mymethod::mutation_best(core::pose::Pose bestpose, core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand3, rand4, rand_fra1, rand_fra2, rand_fra3;;
  core::pose::Pose pose;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1;while(rand_fra3 ==rand_fra2 || rand_fra3 == rand_fra1);
   
  pose = bestpose;
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    pose.set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra2 + k, mypose[rand2].phi(rand_fra2 + k));
    pose.set_psi(rand_fra2 + k, mypose[rand2].psi(rand_fra2 + k));
  }
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra3 + k, mypose[rand3].phi(rand_fra3 + k));
    pose.set_psi(rand_fra3 + k, mypose[rand3].psi(rand_fra3 + k));
  }  
  posemutation = pose;
}

void Mymethod::convex_mutation3(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand_fra1, rand_fra2;
  core::pose::Pose pose;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);

  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   
  pose = mypose[target_id];
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    pose.set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra2 + k, mypose[rand2].phi(rand_fra2 + k));
    pose.set_psi(rand_fra2 + k, mypose[rand2].psi(rand_fra2 + k));
  } 
  posemutation = pose;
}

void Mymethod::convex_mutation4(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, randElow,rand_fra1, rand_fra2;
  core::pose::Pose pose;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  if(int(ElowerthanCtarget.size()) == 0){
    do randElow = rand()%pop_size;while(randElow == target_id || randElow == rand1 || randElow == rand2);
    pose = mypose[randElow];
  }else{
  
    randElow = rand()%int(ElowerthanCtarget.size());
    pose = ElowerthanCtarget[randElow];
  }
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   

  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    pose.set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose.set_phi(rand_fra2 + k, mypose[rand2].phi(rand_fra2 + k));
    pose.set_psi(rand_fra2 + k, mypose[rand2].psi(rand_fra2 + k));
  } 
  posemutation = pose;
}

void Mymethod::convex_mutation5(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, randElow,rand_fra1, rand_fra2;
  core::pose::Pose poseSL;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  if(int(ElowerthanCtarget.size()) == 0){
    do randElow = rand()%pop_size;while(randElow == target_id || randElow == rand1 );
    poseSL = mypose[randElow];
  }
  else{
    randElow = rand()%int(ElowerthanCtarget.size());
    poseSL = ElowerthanCtarget[randElow];
  }
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   

  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    mypose[target_id].set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra2 + k, poseSL.phi(rand_fra2 + k));
    mypose[target_id].set_psi(rand_fra2 + k, poseSL.psi(rand_fra2 + k));
  } 
  posemutation = mypose[target_id];
}

void Mymethod::convex_mutation6(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, randElow,rand_fra1, rand_fra2;
  core::pose::Pose poseSL;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  
  if(int(ElowerthanCtarget.size()) == 0){
    do randElow = rand()%pop_size;while(randElow == target_id || randElow == rand1 || randElow == rand2);
    poseSL = mypose[randElow];
  }else{
  
    randElow = rand()%int(ElowerthanCtarget.size());
    poseSL = ElowerthanCtarget[randElow];
  }
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   

  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra1 + k, mypose[rand2].phi(rand_fra1 + k));
    mypose[rand1].set_psi(rand_fra1 + k, mypose[rand2].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra2 + k, poseSL.phi(rand_fra2 + k));
    mypose[rand1].set_psi(rand_fra2 + k, poseSL.psi(rand_fra2 + k));
  } 
  posemutation = mypose[rand1];
}

void Mymethod::convex_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length){
  int cr_pos;
  cr_pos = rand()%(seq_length-3)+1;
    
  for(int k = 0 ; k<3 ; k++){
    posemutation.set_phi(cr_pos + k, target.phi(cr_pos + k));
    posemutation.set_psi(cr_pos + k, target.psi(cr_pos + k));
  }
}
 
void Mymethod::SaDE_mutation1(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id){
  int rand1, rand2, rand3, rand_fra1, rand_fra2;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  pose3 = mypose[rand3];
  for(int k = 0 ; k<3 ; k++){
    pose3.set_phi(rand_fra1 + k,pose1.phi(rand_fra1 + k));
    pose3.set_psi(rand_fra1 + k,pose1.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose3.set_phi(rand_fra2 + k,pose2.phi(rand_fra2 + k));
    pose3.set_psi(rand_fra2 + k,pose2.psi(rand_fra2 + k));
  }
  posemutation = pose3; 
}

void Mymethod::SaDE_mutation2(core::pose::Pose mypose[], core::pose::Pose& posemutation, core::pose::Pose poselowE, int pop_size, int seq_length, int target_id){
  int rand1, rand2,rand_fra1, rand_fra2, rand_fra3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1;while(rand_fra3 ==rand_fra2 || rand_fra3 == rand_fra1);
   

  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra1 + k, mypose[rand1].phi(rand_fra1 + k));
    mypose[target_id].set_psi(rand_fra1 + k, mypose[rand1].psi(rand_fra1 + k));
  }
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra2 + k, mypose[rand2].phi(rand_fra2 + k));
    mypose[target_id].set_psi(rand_fra2 + k, mypose[rand2].psi(rand_fra2 + k));
  }  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra3 + k, poselowE.phi(rand_fra3 + k));
    mypose[target_id].set_psi(rand_fra3 + k, poselowE.psi(rand_fra3 + k));
  }   
  posemutation = mypose[target_id];  
}

void Mymethod::SaDE_mutation3(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id){
  int rand1, rand2, rand3, rand4, rand_fra1, rand_fra2, rand_fra3;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2);
  do rand4 = rand()%pop_size; while(rand4 == target_id || rand4 == rand1 || rand4 == rand2 || rand4 == rand3);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1; while(rand_fra3 == rand_fra1 || rand_fra3 == rand_fra2);
   
  pose1 = mypose[rand1];
  pose2 = mypose[rand2];
  pose3 = mypose[rand3];
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra1 + k,pose2.phi(rand_fra1 + k));
    pose1.set_psi(rand_fra1 + k,pose2.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra2 + k,pose3.phi(rand_fra2 + k));
    pose1.set_psi(rand_fra2 + k,pose3.psi(rand_fra2 + k));
  }
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra3 + k,mypose[rand4].phi(rand_fra3 + k));
    pose1.set_psi(rand_fra3 + k,mypose[rand4].psi(rand_fra3 + k));
  } 
  posemutation = pose1; 
  
}

void Mymethod::SaDE_mutation4(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id){
  int rand1, rand2, rand_fra1, rand_fra2;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
   

  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra1 + k,mypose[rand1].phi(rand_fra1 + k));
    mypose[target_id].set_psi(rand_fra1 + k,mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra2 + k,mypose[rand2].phi(rand_fra2 + k));
    mypose[target_id].set_psi(rand_fra2 + k,mypose[rand2].psi(rand_fra2 + k));
  }
  posemutation = mypose[target_id]; 
  
}


void Mymethod::UMDE_mutation1(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id){
  int randbetter, rand1, rand2, rand_fra1, rand_fra2, range;
  range = pop_size * 0.05;
  randbetter = rand()%range; 
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == sort_index[randbetter]);
  do rand2 = rand()%pop_size; while(rand2 == target_id ||  rand2 == rand1 || rand2 == sort_index[randbetter]);

  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  
  for(int k = 0 ; k<3 ; k++){
    sortpose[randbetter].set_phi(rand_fra1 + k,mypose[rand1].phi(rand_fra1 + k));
    sortpose[randbetter].set_psi(rand_fra1 + k,mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    sortpose[randbetter].set_phi(rand_fra2 + k,mypose[rand2].phi(rand_fra2 + k));
    sortpose[randbetter].set_psi(rand_fra2 + k,mypose[rand2].psi(rand_fra2 + k));
  } 
  posemutation = sortpose[randbetter];
}

void Mymethod::UMDE_mutation2(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id){
  int randbetter, rand1, rand2, rand3, rand_fra1, rand_fra2, rand_fra3, range;
  range = pop_size * 0.05;
  randbetter = rand()%range; 
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == sort_index[randbetter]);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1 || rand2  == sort_index[randbetter]);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2 || rand3 == sort_index[randbetter]);

  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1; while(rand_fra3 == rand_fra1 || rand_fra3 == rand_fra2);
  
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra1 + k,sortpose[randbetter].phi(rand_fra1 + k));
    mypose[rand1].set_psi(rand_fra1 + k,sortpose[randbetter].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra2 + k,mypose[rand2].phi(rand_fra2 + k));
    mypose[rand1].set_psi(rand_fra2 + k,mypose[rand2].psi(rand_fra2 + k));
  } 
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra3 + k,mypose[rand3].phi(rand_fra3 + k));
    mypose[rand1].set_psi(rand_fra3 + k,mypose[rand3].psi(rand_fra3 + k));
  }  
  posemutation = mypose[rand1];
}

void Mymethod::UMDE_mutation3(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id){
  int randbetter, rand1, rand_fra1, rand_fra2, range;
  range = pop_size * 0.05;
  randbetter = rand()%range; 
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == sort_index[randbetter]);

  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra1 + k,sortpose[randbetter].phi(rand_fra1 + k));
    mypose[target_id].set_psi(rand_fra1 + k,sortpose[randbetter].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra2 + k,mypose[rand1].phi(rand_fra2 + k));
    mypose[target_id].set_psi(rand_fra2 + k,mypose[rand1].psi(rand_fra2 + k));
  } 
  posemutation = mypose[target_id];
}

void Mymethod::UMDE_mutation4(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index){
  int rand1, rand2, rand_fra1, rand_fra2;
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == index);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1 || rand2 == index);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  
  for(int k = 0 ; k<3 ; k++){
    bestpose.set_phi(rand_fra1 + k,mypose[rand1].phi(rand_fra1 + k));
    bestpose.set_psi(rand_fra1 + k,mypose[rand1].psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    bestpose.set_phi(rand_fra2 + k,mypose[rand2].phi(rand_fra2 + k));
    bestpose.set_psi(rand_fra2 + k,mypose[rand2].psi(rand_fra2 + k));
  } 
  posemutation = bestpose;
}
void Mymethod::UMDE_mutation5(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index){
  int rand1, rand2, rand3, rand_fra1, rand_fra2, rand_fra3;
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == index);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1 || rand2 == index);
  do rand3 = rand()%pop_size; while(rand3 == target_id || rand3 == rand1 || rand3 == rand2 || rand3 == index);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  do rand_fra3 = rand()%(seq_length-3)+1; while(rand_fra3 == rand_fra1 || rand_fra3 == rand_fra2);
  
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra1 + k,bestpose.phi(rand_fra1 + k));
    mypose[rand1].set_psi(rand_fra1 + k,bestpose.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra2 + k,mypose[rand2].phi(rand_fra2 + k));
    mypose[rand1].set_psi(rand_fra2 + k,mypose[rand2].psi(rand_fra2 + k));
  } 
  for(int k = 0 ; k<3 ; k++){
    mypose[rand1].set_phi(rand_fra3 + k,mypose[rand3].phi(rand_fra3 + k));
    mypose[rand1].set_psi(rand_fra3 + k,mypose[rand3].psi(rand_fra3 + k));
  }   
  posemutation = mypose[rand1];
}

void Mymethod::UMDE_mutation6(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index){
  int rand1, rand_fra1, rand_fra2;
  do rand1 = rand()%pop_size; while(rand1 == target_id || rand1 == index);
  
  rand_fra1 = rand()%(seq_length-3)+1;
  do rand_fra2 = rand()%(seq_length-3)+1; while(rand_fra2 == rand_fra1);
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra1 + k,bestpose.phi(rand_fra1 + k));
    mypose[target_id].set_psi(rand_fra1 + k,bestpose.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    mypose[target_id].set_phi(rand_fra2 + k,mypose[rand1].phi(rand_fra2 + k));
    mypose[target_id].set_psi(rand_fra2 + k,mypose[rand1].psi(rand_fra2 + k));
  } 
  posemutation = mypose[target_id];
}


void Mymethod::zxg_mutation(core::pose::Pose mypose[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id){
  int rand1, rand2, rand_fra1, rand_fra2, fra_pos1, fra_pos2;
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
  
  
  fra_pos1 = rand()%(seq_length - 4) + 1;
  do fra_pos2 = rand()%(seq_length - 4) + 1; while(fra_pos2 == fra_pos1);
  
  rand_fra1 = rand()%17 + 3;
  do rand_fra2 = rand()%17 + 3; while(rand_fra2 == rand_fra1);
  
  posemutation = mypose[target_id];
  
  for(int k = 0 ; k < rand_fra1; k++){
    if(fra_pos1 + k <= seq_length){
      posemutation.set_phi(fra_pos1 + k,mypose[rand1].phi(fra_pos1 + k));
      posemutation.set_psi(fra_pos1 + k,mypose[rand1].phi(fra_pos1 + k));
    }
  }
  
  for(int k = 0 ; k < rand_fra2; k++){
    if(fra_pos2 + k <= seq_length){
      posemutation.set_phi(fra_pos2 + k,mypose[rand2].phi(fra_pos2 + k));
      posemutation.set_psi(fra_pos2 + k,mypose[rand2].phi(fra_pos2 + k));
    }
  }
}

void Mymethod::MutationandCrossover_RandomSelectedFrag(core::pose::Pose mypose[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id){
  int rand1, rand2, rand_fra1, rand_fra2, fra_pos1, fra_pos2, cr_pos;
  double probabilityofCrossover;
  do rand1 = rand()%pop_size; while(rand1 == target_id);
//  do rand2 = rand()%pop_size; while(rand2 == target_id || rand2 == rand1);
    
  rand_fra1 = rand()%7 + 3;
  
  fra_pos1 = rand()%(seq_length - rand_fra1) + 1;
//  do fra_pos2 = rand()%(seq_length - 4) + 1; while(fra_pos2 == fra_pos1);
  

//  do rand_fra2 = rand()%17 + 3; while(rand_fra2 == rand_fra1);
  
  posemutation = mypose[target_id];
  
  for(int k = 0 ; k < rand_fra1; k++){
    if(fra_pos1 + k <= seq_length){
      posemutation.set_phi(fra_pos1 + k,mypose[rand1].phi(fra_pos1 + k));
      posemutation.set_psi(fra_pos1 + k,mypose[rand1].phi(fra_pos1 + k));
    }
  }
  
 /* for(int k = 0 ; k < rand_fra2; k++){
    if(fra_pos2 + k <= seq_length){
      posemutation.set_phi(fra_pos2 + k,mypose[rand2].phi(fra_pos2 + k));
      posemutation.set_psi(fra_pos2 + k,mypose[rand2].phi(fra_pos2 + k));
    }
  }*/
 
 //////////////////////////////////////////////////////////////////////////crossover
  probabilityofCrossover = rand()%100/(double)101;
  if(probabilityofCrossover > 0.5){
	cr_pos = rand()%(seq_length - rand_fra1) + 1;
	for(int k = 0 ; k < rand_fra1; k++){
		if(cr_pos + k <= seq_length){
			posemutation.set_phi(cr_pos + k,mypose[target_id].phi(cr_pos + k));
			posemutation.set_psi(cr_pos + k,mypose[target_id].phi(cr_pos + k));
		}
	} 
  }
}

void Mymethod::mutation_strategy1_Metacontact(core::pose::Pose posevector[],core::pose::Pose lowcontact_pose2,core::pose::Pose lowcontact_pose3,core::pose::Pose lowcontact_pose4,
				  core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand_fra1, rand_fra2, randpose;
  core::pose::Pose pose1, pose2;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  
  randpose = rand()%3 + 1;
  
  
  rand_fra1 = rand()%(seq_length-9)+1;
  do rand_fra2 = rand()%(seq_length-9)+1; while(rand_fra2 == rand_fra1);
   
  posemutation = posevector[target_id];
  pose1 = posevector[rand1];
  
  switch (randpose){
    case 1 :
      pose2 = lowcontact_pose2;
      break;
    case 2 :
      pose2 = lowcontact_pose3;
      break;
    case 3 :
      pose2 = lowcontact_pose4;
      break;
    default :
      pose2 = lowcontact_pose2;     
  }
  
  
  for(int k = 0 ; k<9 ; k++){
    posemutation.set_phi(rand_fra1 + k, pose1.phi(rand_fra1 + k));
    posemutation.set_psi(rand_fra1 + k, pose1.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<9 ; k++){
    posemutation.set_phi(rand_fra2 + k, pose2.phi(rand_fra2 + k));
    posemutation.set_psi(rand_fra2 + k, pose2.psi(rand_fra2 + k));
  }  
}

void Mymethod::mutation_strategy1(core::pose::Pose posevector1[],core::pose::Pose posevector2[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],
				  core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id)
{
  int rand1, rand2, rand3, rand_fra1, rand_fra2, rand_pop1, rand_pop2;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand1 = rand()%pop_size; while(rand1 == target_id);
  
  rand_pop1 = rand()%3 + 1;
  do rand_pop2 = rand()%3 + 1; while(rand_pop2 == rand_pop1);
  
  rand2 = rand()%pop_size;
  rand3 = rand()%pop_size;
  
  rand_fra1 = rand()%(seq_length-9)+1;
  do rand_fra2 = rand()%(seq_length-9)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = posevector1[rand1];
  
  switch (rand_pop1){
    case 1 :
      pose2 = posevector2[rand2];
      break;
    case 2 :
      pose2 = posevector3[rand2];
      break;
    case 3 :
      pose2 = posevector4[rand2];
      break;
    default :
      pose2 = posevector2[rand2];     
  }
  
   switch (rand_pop2){
    case 1 :
      pose3 = posevector2[rand3];
      break;
    case 2 :
      pose3 = posevector3[rand3];
      break;
    case 3 :
      pose3 = posevector4[rand3];
      break;
    default :
      pose3 = posevector2[rand3];     
  }
  
  for(int k = 0 ; k<9 ; k++){
    pose1.set_phi(rand_fra1 + k,pose2.phi(rand_fra1 + k));
    pose1.set_psi(rand_fra1 + k,pose2.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra2 + k,pose3.phi(rand_fra2 + k));
    pose1.set_psi(rand_fra2 + k,pose3.psi(rand_fra2 + k));
  }  
  posemutation = pose1;
}

void Mymethod::mutation_strategy2(core::pose::Pose posevector1[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],core::pose::Pose bestpose,
				  core::pose::Pose& posemutation,int pop_size,int seq_length)
{
  int rand2, rand3, rand_fra1, rand_fra2, rand_pop1, rand_pop2;
  core::pose::Pose pose1, pose2, pose3;
  
  do rand_pop1 = rand()%4; while(rand_pop1 == 1);
  do rand_pop2 = rand()%4; while(rand_pop2 == rand_pop1 || rand_pop2 == 1);
  
  rand2 = rand()%pop_size;
  rand3 = rand()%pop_size;
  
  rand_fra1 = rand()%(seq_length-9)+1;
  do rand_fra2 = rand()%(seq_length-9)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = bestpose;
  
  switch (rand_pop1){
    case 0 :
      pose2 = posevector1[rand2];
      break;
    case 2 :
      pose2 = posevector3[rand2];
      break;
    case 3 :
      pose2 = posevector4[rand2];
      break;
    default :
      pose2 = posevector1[rand2];     
  }
  
   switch (rand_pop2){
    case 0 :
      pose3 = posevector1[rand3];
      break;
    case 2 :
      pose3 = posevector3[rand3];
      break;
    case 3 :
      pose3 = posevector4[rand3];
      break;
    default :
      pose3 = posevector3[rand3];     
  }
  
  for(int k = 0 ; k<9 ; k++){
    pose1.set_phi(rand_fra1 + k,pose2.phi(rand_fra1 + k));
    pose1.set_psi(rand_fra1 + k,pose2.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra2 + k,pose3.phi(rand_fra2 + k));
    pose1.set_psi(rand_fra2 + k,pose3.psi(rand_fra2 + k));
  }  
  posemutation = pose1;

}

void Mymethod::mutation_strategy3(core::pose::Pose posevector1[],core::pose::Pose posevector2[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],
				  core::pose::Pose& posemutation,int pop_size,int seq_length)
{
  int rand1, rand2, rand3, rand_fra1, rand_fra2, rand_pop1, rand_pop2;
  core::pose::Pose pose1, pose2, pose3;
  
  rand1 = rand()%(pop_size/2);
  
  do rand_pop1 = rand()%4; while(rand_pop1 == 2);
  do rand_pop2 = rand()%4; while(rand_pop2 == rand_pop1 || rand_pop1 == 2);
  
  rand2 = rand()%pop_size;
  rand3 = rand()%pop_size;
  
  rand_fra1 = rand()%(seq_length-9)+1;
  do rand_fra2 = rand()%(seq_length-9)+1; while(rand_fra2 == rand_fra1);
   
  pose1 = posevector3[rand1];
  
  switch (rand_pop1){
    case 0 :
      pose2 = posevector1[rand2];
      break;
    case 1 :
      pose2 = posevector2[rand2];
      break;
    case 3 :
      pose2 = posevector4[rand2];
      break;
    default :
      pose2 = posevector1[rand2];     
  }
  
   switch (rand_pop2){
    case 0 :
      pose3 = posevector2[rand3];
      break;
    case 1 :
      pose3 = posevector3[rand3];
      break;
    case 3 :
      pose3 = posevector4[rand3];
      break;
    default :
      pose3 = posevector2[rand3];     
  }
  
  for(int k = 0 ; k<9 ; k++){
    pose1.set_phi(rand_fra1 + k,pose2.phi(rand_fra1 + k));
    pose1.set_psi(rand_fra1 + k,pose2.psi(rand_fra1 + k));
  }
  
  for(int k = 0 ; k<3 ; k++){
    pose1.set_phi(rand_fra2 + k,pose3.phi(rand_fra2 + k));
    pose1.set_psi(rand_fra2 + k,pose3.psi(rand_fra2 + k));
  }  
  posemutation = pose1;
}


float Mymethod::DPerror(std::vector<int>& pd_position_first, std::vector<int>& pd_position_second, std::vector<int>& pd_peak_number, std::vector<float>& pd_peak_distance, std::vector<float>& decoy_ditance){

  float DPscore = 0.0;
  float D = 0.01;
  float F = 0.01;
  for(int i = 0; i < int(pd_position_first.size()); i++){
  
    D = fabs(decoy_ditance[i] - pd_peak_distance[i]);
    F = pd_peak_number[i];
    DPscore =DPscore + F/D;
  }
  return DPscore;
}

double Mymethod::matchingScore_of_Saxs(std::string path){
	 
	 std::string line;
	 double q = 0.0;
	 double exp_intensity = 0.0;
	 double error = 0.0;
	 double model_intensity = 0.0;
	 double SAXS_Score = 0.0;
	 std::ifstream openfile(path.c_str());
	 while(getline(openfile, line)){
		if(line.substr(0,1) != "#"){
			std::istringstream is(line);
			is>> q>> exp_intensity>> error>> model_intensity;		  
			SAXS_Score = SAXS_Score + ((exp_intensity - model_intensity) * (exp_intensity - model_intensity)) * q;			
		}
	 }
	 openfile.close();
	 openfile.clear();
	 return SAXS_Score;
}
double Mymethod::get_chi_fromsaxs(std::string path){
	 std::string line;
	 std::string chi_string;
	 double chi = 0.0;
	 int flag = 0;
	 
	 std::ifstream openfile(path.c_str());
	 while(getline(openfile, line) && flag < 2){
		if(flag == 1){
			std::istringstream is(line);
			for(int m =0; m < 11; m++)
				is>>chi_string;
		}
		chi = atof(chi_string.c_str());
		flag++;		
	 }
	 openfile.close();
	 openfile.clear();
	 return chi;
}

double Mymethod::contactScore(core::pose::Pose pose, vector<contactofTop> vectorContact, int seq_length){

	double Econt = 0.0;
	double residues_distance = 0.0;
	double Uij = 0.0;
	for(int i =0; i < int(vectorContact.size()); i++){
		residues_distance = core::scoring::residues_distance_CB(pose, vectorContact[i].first_residue, vectorContact[i].second_residue);
		Uij = vectorContact[i].contactofij;
		if(residues_distance < 8)
			Econt = Econt - Uij;
		else if(residues_distance < 10)
			Econt = Econt - 0.5*Uij*(1 - sin(((residues_distance - 9)/2)*3.14));
			else if(residues_distance < 80)
				Econt = Econt - 0.5*Uij*(1 + sin(((residues_distance - 45)/70)*3.14));
				else
					Econt = Econt + Uij;
	}
	return Econt;
}






