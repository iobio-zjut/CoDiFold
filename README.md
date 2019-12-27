# CoDiFold
### a de novo protein structure prediction by coupling contact with distance profile

**Developer:**   
                Chunxiang Peng  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: pengcx@zjut.edu.cn  
		
**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure CoDiFold:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/ 
and extract it to ``"~/"`` directory.

- Copy and paste source code of ``"ClassicAbinitio.hh"``, ``"ClassicAbinitio.cc"``, ``"ClassicAbinitio.fwd.hh"``,  ``"mymethod.hh"``,  ``"mymethod.cc"``, ``"readfile.hh"`` and ``"readfile.cc"`` from ``"src/"`` folder in CoDiFold package to ``"~/Rosetta/main/source/src/protocols/abinitio/"`` folder in Rosetta. Copy and paste source code of ``"rms_util.hh"``, ``"rms_util.cc"`` and ``"rms_util.tmpl.hh"`` from ``"src/"`` folder in CoDiFold package to ``"~/Rosetta/main/source/src/core/scoring/"`` folder in Rosetta. Copy and paste configuration file ``"protocols_b_6.src.settings"`` from ``"src/"`` folder in CoDiFold package to ``"~/Rosetta/main/source/src/"`` folder in Rosetta.

- Compile Rosetta source code using the following commands:

```
 $ cd ~/Rosetta/main/source/
 $ ./scons.py -j<NumOfJobs> mode=release bin
```

- If you want to recompile CoDiFold source code, use the following commands:

```
 $ cd ~/CoDiFoldFold/
 $ g++ -o bin/CoDiFold src/CoDiFold.cpp
```
## 2. INPUT
CoDiFold requires six files to generate models:

	-f	seq.fasta		: protein sequence fasta file
	-c1	contact1.txt		: contact map file from ResPRE
	-c2	contact2.txt		: contact map file from Raptor-X
	-dp	dp.txt			: distance profile file
	-frag3	3mer_fragment_library	: fragment library with fragment lenth 3
	-frag9	9mer_fragment_library	: fragment library with fragment lenth 9

## 3. OUTPUT
Output file of CoDiFold is stored in current folder, including one predicted model (lowestEtotal.pdb).

	lowestEtotal.pdb	the lowest energy model in whole process

## 4. EXAMPLE
Please follow the below steps to run CoDiFold:

- Go to the ``"example/"`` folder of CoDiFold.
  
- Run CoDiFold with the following command:
  
```
 $ cd ~/CoDiFold/
 $ ./../bin/CoDiFold -f ./input_files/seq.fasta -c1 input_files/contact1.txt -c2 input_files/contact2.txt -dp input_files/dp.txt -frag3 input_files/3mer_fragment_library -frag9 input_files/9mer_fragment_library
```

- Final model is generated in current folder.

## 5. DISCLAIMER
The executable software and the source code of CoDiFold is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
