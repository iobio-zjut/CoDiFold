//
//  extract_dp.hpp
//  contact
//
//  Created by 彭春祥 on 2018/6/16.
//  Copyright © 2018年 彭春祥. All rights reserved.
//


#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>

#include <utility/vector1.hh>

/////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
#ifndef readfile_hh
#define readfile_hh
#endif /* readfile_hh */

using namespace std;

struct contactofTop{
	int first_residue;
	int second_residue;
	double contactofij;
};

bool comp(const contactofTop &a, const contactofTop &b);

/////////////////////////////////////////////////////////////////////////////

void readcontact(vector<contactofTop>& vectorContact, string contactName);
void extract_distance_profile(vector<int>& pd_position_first,vector<int>& pd_position_second,vector<int>& pd_totalNumber,vector<int>& pd_peak_number,vector <float>& pd_peak_distance);
void extract_distance_number(vector<int>& distance_number375,vector<int>& distance_number425,vector<int>& distance_number475,vector<int>& distance_number525,vector<int>& distance_number575
			      ,vector<int>& distance_number625,vector<int>& distance_number675,vector<int>& distance_number725,vector<int>& distance_number775,vector<int>& distance_number825
			      ,vector<int>& distance_number875,vector<int>& distance_number925);
void fun_extract_contact1(double **arr, int seq_length);
void fun_extract_contact2(double **arr, int seq_length);
double **new2DDoubleArr(int row, int col);
void release2DArr(int n, double **Arr);
                       		   	
