//
//  extract_file.cpp
//  contact
//
//  Created by 彭春祥 on 2018/6/16.
//  Copyright © 2018年 彭春祥. All rights reserved.
//

#include "readfile.hh"


// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>
#include <core/id/types.hh>
#include <core/chemical/automorphism.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/TMscore.hh>

// Utility headers
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

///////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <string.h>	
#include <set>
#include <string>
#include <algorithm>
#include <stdlib.h>
//////////////////////////////////////////////////////////

using namespace std;



void readcontact(vector<contactofTop>& vectorContact, string contactName){

//	int NumberofContacts = (int)0.2 * seq_length;
	string line;
	ifstream myfile1(contactName.c_str());
	while(getline (myfile1, line)){
		contactofTop onecontact;
		int rsd_start,rsd_end;
		double contact = 0.0;
		istringstream is(line);
		is >> rsd_start >> rsd_end >> contact;
		if( fabs(rsd_end - rsd_start) >= 6){
			onecontact.first_residue = rsd_start;
			onecontact.second_residue = rsd_end;
			onecontact.contactofij = contact;
			vectorContact.push_back(onecontact);
		}
	}
	
	sort(vectorContact.begin(), vectorContact.end(), comp);
	myfile1.close();
	myfile1.clear();
}

void fun_extract_contact1(double **arr, int seq_length)//  read contact
{
    int NumberofContacts = 0.5 * seq_length;
    int count = 0;
    string line;
    ifstream myfile1("input_files/contact1.txt");  
    
    if(myfile1.is_open())
      std::cout<<"contact1 open success!!!"<<std::endl;
    else
      std::cout<<"contact1 open failed!!!!"<<std::endl;
    
    while(getline (myfile1, line)) //  && count < NumberofContacts
      {
	istringstream is(line);
	int rsd_start,rsd_end;
	double contact_p;
	is>>rsd_start>>rsd_end>>contact_p;
	if( fabs(rsd_end - rsd_start) >= 6 && contact_p >= 0.5){
	  
	  arr[rsd_start-1][rsd_end-1] = contact_p;
	  count++;
	}
      }
    myfile1.close();
    myfile1.clear();	  
	  
}

void fun_extract_contact2(double **arr, int seq_length)//  read contact	
{
  int NumberofContacts = 0.5 * seq_length;
  int count = 0;
  string line;
  ifstream myfile2("input_files/contact2.txt");
  
  if(myfile2.is_open())
	std::cout<<"contact2 open success!!!"<<std::endl;
  else
	std::cout<<"contact2 open failed!!!!"<<std::endl;
  
  while(getline (myfile2, line))   //  && count < NumberofContacts
    {
       if ( !((line[0] >= 'A' && line[0] <= 'Z') || (line[0] >= 'a' && line[0] <= 'z')) ) {
	  istringstream is(line);
	  int rsd_start,rsd_end,str0,str8;
	  double contact_p;
	  is>>rsd_start>>rsd_end>>str0>>str8>>contact_p;
	  if( fabs(rsd_end - rsd_start) >= 6 && contact_p >= 0.5){
	  arr[rsd_start-1][rsd_end-1] = contact_p;
	  count++;
	}
       }    
    }  
  myfile2.close();
  myfile2.clear();                                          
}

void extract_distance_number(vector<int>& distance_number375,vector<int>& distance_number425,vector<int>& distance_number475,vector<int>& distance_number525,vector<int>& distance_number575
			      ,vector<int>& distance_number625,vector<int>& distance_number675,vector<int>& distance_number725,vector<int>& distance_number775,vector<int>& distance_number825
			      ,vector<int>& distance_number875,vector<int>& distance_number925){
    ifstream myfile("input_files/dp.txt");
    if (!myfile.is_open()){
        cout << "Unable to open myfile";
    }
    string temp;
    string *arr_dp;
    arr_dp = new string[26];
    while (getline(myfile,temp)) {
      if(temp.size() > 5){
        istringstream is(temp);
        
        is>>arr_dp[0]>>arr_dp[1]>>arr_dp[2]>>arr_dp[3]>>arr_dp[4]>>arr_dp[5]>>arr_dp[6]>>arr_dp[7]>>arr_dp[8]>>arr_dp[9]>>arr_dp[10]>>arr_dp[11]>>arr_dp[12]>>arr_dp[13]>>arr_dp[14]>>arr_dp[15]>>arr_dp[16]>>arr_dp[17]>>arr_dp[18]>>arr_dp[19]>>arr_dp[20]>>arr_dp[21]>>arr_dp[22]>>arr_dp[23]>>arr_dp[24]>>arr_dp[25];
	distance_number375.push_back(atoi(arr_dp[13].c_str()));
	distance_number425.push_back(atoi(arr_dp[14].c_str()));
	distance_number475.push_back(atoi(arr_dp[15].c_str()));
	distance_number525.push_back(atoi(arr_dp[16].c_str()));
	distance_number575.push_back(atoi(arr_dp[17].c_str()));
	distance_number625.push_back(atoi(arr_dp[18].c_str()));
	distance_number675.push_back(atoi(arr_dp[19].c_str()));
	distance_number725.push_back(atoi(arr_dp[20].c_str()));
	distance_number775.push_back(atoi(arr_dp[21].c_str()));
	distance_number825.push_back(atoi(arr_dp[22].c_str()));
	distance_number875.push_back(atoi(arr_dp[23].c_str()));
	distance_number925.push_back(atoi(arr_dp[24].c_str()));
      }
    } 
    delete []arr_dp;
    myfile.clear();
    myfile.close();
}
void extract_distance_profile(vector<int>& pd_position_first,vector<int>& pd_position_second,vector<int>& pd_totalNumber,vector<int>& pd_peak_number,vector <float>& pd_peak_distance){                  //     read distance                   
    ifstream myfile("input_files/dp.txt");
    if (!myfile.is_open()){
        cout << "Unable to open myfile";
  //      system("pause");
     //   exit(1);
    }
    
       string temp;

    while (getline(myfile,temp)) {
      if(temp.size() > 5){
        istringstream is(temp);
        string arr_dp[26];
        is>>arr_dp[0]>>arr_dp[1]>>arr_dp[2]>>arr_dp[3]>>arr_dp[4]>>arr_dp[5]>>arr_dp[6]>>arr_dp[7]>>arr_dp[8]>>arr_dp[9]>>arr_dp[10]>>arr_dp[11]>>arr_dp[12]>>arr_dp[13]>>arr_dp[14]>>arr_dp[15]>>arr_dp[16]>>arr_dp[17]>>arr_dp[18]>>arr_dp[19]>>arr_dp[20]>>arr_dp[21]>>arr_dp[22]>>arr_dp[23]>>arr_dp[24]>>arr_dp[25];
        pd_position_first.push_back(atoi(arr_dp[0].c_str()));
        pd_position_second.push_back(atoi(arr_dp[1].c_str()));
        pd_totalNumber.push_back(atoi(arr_dp[4].c_str()));
        pd_peak_number.push_back(atoi(arr_dp[24].c_str()));
        pd_peak_distance.push_back(atof(arr_dp[25].c_str()));
      }
    }
    myfile.clear();
    myfile.close();
}

double **new2DDoubleArr(int row, int col){

  double **ans = new double*[row];
  for(int i = 0; i< row; i++){
  
    ans[i] = new double[col];
  }
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++)  
      ans[i][j] = 0.0;
  }
  return ans;
}

void release2DArr(int n, double **Arr){

  for(int i = 0; i < n; i++){
  
    delete [] Arr[i];
  }
  delete [] Arr;
}

bool comp(const contactofTop &a, const contactofTop &b){
    return a.contactofij > b.contactofij;
}

