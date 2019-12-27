

#include <string>
#include<stdio.h>
#include <vector>
#include <string>
#include<math.h>

// project
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include "readfile.hh"


class Mymethod{
  

public:
  
  std::vector<double> coefficient_contact;
  float sum_energy_distanceprofile;
  int both_exist_number;
  int pop_size;
  int belta;
  float CR;
  

  void computed_coefficient_contact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double **contact1,double **contact2);
  void energy_distanceprofile(std::vector<float>& coefficient_contact,std::vector<float>& decoy_ditance,std::vector <float>& pd_peak_distance);
  void computed_bothexist_dpandcontact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double contact1[][108],double contact2[][108]);
  void computed_distanceprofile(std::vector<float>& decoy_ditance,std::vector <float>& pd_peak_distance);
  void show_dppair_contact(std::vector<int>& pd_position_first,std::vector<int>& pd_position_second,double contact1[][108],double contact2[][108]);
  void dp_initialization(core::pose::Pose& pose,core::pose::Pose mypose[] );
  void dp_mutation(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void dp_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length);
  bool average_distance_error_discriminate(const core::pose::Pose& pose_trila,const core::pose::Pose& pose_target,std::vector<int>& pd_position_first,
			    std::vector<int>& pd_position_second,std::vector<float>& pd_peak_distance,std::vector<float>& decoy_ditance_trial,std::vector<float>& decoy_ditance_target);
  
  void energy_distanceprofile2(std::vector<float>& decoy_ditance,std::vector<float>& coefficient_contact,float& sumE,std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925);
  void energy_distanceprofile3(std::vector<float>& decoy_ditance,std::vector<double>& coefficient_contact,double& sumE, std::vector<int>& pd_totalNumber,std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925);
  void energy3(std::vector<float>& decoy_ditance,std::vector<int>& pd_totalNumber,float& sumE,std::vector<float>& delta_E,std::vector<int>& distance_number375,std::vector<int>& distance_number425,std::vector<int>& distance_number475,std::vector<int>& distance_number525,std::vector<int>& distance_number575
			      ,std::vector<int>& distance_number625,std::vector<int>& distance_number675,std::vector<int>& distance_number725,std::vector<int>& distance_number775,std::vector<int>& distance_number825
			      ,std::vector<int>& distance_number875,std::vector<int>& distance_number925);
  bool dp_contact_select(std::vector<float>& dp_E_trial, std::vector<float>& dp_E_target, float sumE_trial, float sumE_target);
  void loop_crossover(core::pose::Pose mypose[], core::pose::Pose& posemutation, std::vector<int>& loop_pos_start, std::vector<int>& loop_pos_end, int pop_size, int target_id);
  void set_default();
  
  
  ////////convex
  double find_the_closest(double trial[], double decoy[][600], int seq_length, int kk);
  double find_the_closest_zxg(double *trial, double **decoy, int seq_length, int kk);
  void convex_mutation1(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_mutation2(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void mutation_best(core::pose::Pose bestpose, core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_mutation3(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_mutation4(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_mutation5(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_mutation6(core::pose::Pose mypose[], std::vector<core::pose::Pose>& ElowerthanCtarget, core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void convex_crossover(core::pose::Pose& target, core::pose::Pose& posemutation, int seq_length);
  
  //////////sade
  void SaDE_mutation1(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void SaDE_mutation2(core::pose::Pose mypose[], core::pose::Pose& posemutation, core::pose::Pose poselowE, int pop_size, int seq_length, int target_id);
  void SaDE_mutation3(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  void SaDE_mutation4(core::pose::Pose mypose[],core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  
  ///////////UMDE
  void UMDE_mutation1(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id);
  void UMDE_mutation2(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id);
  void UMDE_mutation3(core::pose::Pose mypose[], core::pose::Pose sortpose[], int sort_index[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id);  
  void UMDE_mutation4(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index);
  void UMDE_mutation5(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index);
  void UMDE_mutation6(core::pose::Pose mypose[], core::pose::Pose bestpose, core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id, int index);
  
  /////////////////////zxg
  void zxg_mutation(core::pose::Pose mypose[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id);
  ////////////////
  void MutationandCrossover_RandomSelectedFrag(core::pose::Pose mypose[], core::pose::Pose& posemutation, int pop_size, int seq_length, int target_id);
  
  void mutation_strategy1(core::pose::Pose posevector1[],core::pose::Pose posevector2[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],
			  core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  
  void mutation_strategy2(core::pose::Pose posevector1[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],core::pose::Pose bestpose,
				  core::pose::Pose& posemutation,int pop_size,int seq_length);
  
  void mutation_strategy3(core::pose::Pose posevector1[],core::pose::Pose posevector2[],core::pose::Pose posevector3[],core::pose::Pose posevector4[],
				  core::pose::Pose& posemutation,int pop_size,int seq_length);
  void mutation_strategy1_Metacontact(core::pose::Pose posevector1[],core::pose::Pose lowcontact_pose2,core::pose::Pose lowcontact_pose3,core::pose::Pose lowcontact_pose4,
				  core::pose::Pose& posemutation,int pop_size,int seq_length,int target_id);
  
  float DPerror(std::vector<int>& pd_position_first, std::vector<int>& pd_position_second, std::vector<int>& pd_peak_number, std::vector<float>& pd_peak_distance, std::vector<float>& decoy_ditance);
  ////////////////////////////////////////////////////////////saxs
  double matchingScore_of_Saxs(std::string path);
  double get_chi_fromsaxs(std::string path);
  
  double contactScore(core::pose::Pose pose, vector<contactofTop> vectorContact, int seq_length);
  
};

  



