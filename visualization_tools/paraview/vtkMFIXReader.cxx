/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// Thanks to Phil Nicoletti and Brian Dotson at the National Energy 
// Technology Laboratory who developed this class.
// Please address all comments to Phil Nicoletti (philip.nicoletti@pp.netl.doe.gov)
// or Brian Dotson (brian.dotson@netl.doe.gov)

#include "vtkMFIXReader.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkErrorCode.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkHexahedron.h"
#include "vtkFloatArray.h"
#include <string>
#include "vtkDataArraySelection.h"
#include "vtkWedge.h"

vtkCxxRevisionMacro(vtkMFIXReader, "$Revision$");
vtkStandardNewMacro(vtkMFIXReader);

//----------------------------------------------------------------------------
vtkMFIXReader::vtkMFIXReader()
{
  this->FileName = NULL;
  this->NumberOfCells = 0;
  this->NumberOfPoints = 0;
  this->NumberOfCellFields = 0;
  this->RequestInformationFlag = 0;
  this->MakeMeshFlag = 0;
  
  this->min = NULL;
  this->max = NULL;
  this->veclen = NULL;
  this->cell_data_array = NULL;
  this->spx_timestep_index_table = NULL;
  
  this->P_ref     = 0;
  this->P_scale   = 1;
  this->DIM_IC    = 5;
  this->DIM_BC    = 5;
  this->DIM_C     = 5;
  this->DIM_IS    = 5;
  this->c_e       = 1.0;
  this->c_f       = 0.0;
  this->phi       = 0.0;
  this->phi_w     = 0.0;
  this->nspx_use  = 9;
  this->NScalar   = 0;
  this->bKepsilon = false;  
  this->nRR       = 0;
  this->nTimes    = 0;
  this->ext[0] = '1';
  this->ext[1] = '2';
  this->ext[2] = '3';
  this->ext[3] = '4';
  this->ext[4] = '5';
  this->ext[5] = '6';
  this->ext[6] = '7';
  this->ext[7] = '8';
  this->ext[8] = '9';
  this->ext[9] = 'A';
  this->ext[10] = 'B';
  this->version_number = 0;

  this->cell_data_array = NULL;
   
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->points = vtkPoints::New();
  this->mesh = vtkUnstructuredGrid::New();
  this->aHexahedron = vtkHexahedron::New();
  this->aWedge = vtkWedge::New();
  this->NMAX = vtkIntArray::New();
  this->C = vtkDoubleArray::New();
  this->DX = vtkDoubleArray::New();
  this->DY = vtkDoubleArray::New();
  this->DZ = vtkDoubleArray::New();
  this->tmpI = vtkIntArray::New();
  this->tmpD = vtkDoubleArray::New();   
  this->FLAG = vtkIntArray::New();
  this->variable_names = vtkStringArray::New();
  this->variable_components = vtkIntArray::New();
  this->variableIndexToSPX = vtkIntArray::New();
  this->var_timesteps = vtkIntArray::New();
  this->var_timestep_table = vtkIntArray::New();
  this->spx_to_nvar_table = vtkIntArray::New();
  this->var_to_skip_table = vtkIntArray::New();
  this->SpxFileExists = vtkIntArray::New();

  this->SetNumberOfInputPorts(0);
  
  // Time support:
  this->TimeStep = 0; // By default the file does not have timestep
  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = 0;
  this->NumberOfTimeSteps = 1;
  this->TimeSteps = 0;
  this->CurrentTimeStep = 0;
  this->TimeStepWasReadOnce = 0;
  
}

//----------------------------------------------------------------------------
vtkMFIXReader::~vtkMFIXReader()
{

  if (this->FileName)
    {
    delete [] this->FileName;
    }

  for (int j=0; j<=this->variable_names->GetMaxId(); j++){
     this->cell_data_array[j]->Delete();
  }
        
  this->CellDataArraySelection->Delete();
  this->points->Delete();
  this->mesh->Delete();
  this->aHexahedron->Delete();
  this->aWedge->Delete();
  this->NMAX->Delete();
  this->C->Delete();
  this->DX->Delete();
  this->DY->Delete();
  this->DZ->Delete();
  this->tmpI->Delete();
  this->tmpD->Delete();
  this->FLAG->Delete();
  this->variable_names->Delete();
  this->variable_components->Delete();
  this->variableIndexToSPX->Delete();
  this->var_timesteps->Delete();
  this->var_timestep_table->Delete();
  this->spx_to_nvar_table->Delete();
  this->var_to_skip_table->Delete();
  this->SpxFileExists->Delete();
  
  if(this->cell_data_array){
  	delete [] this->cell_data_array;
  }
  if(this->min){
  	delete [] this->min;
  }
  if(this->max){
  	delete [] this->max;
  }
  if(this->veclen){
  	delete [] this->veclen;
  }
  if(this->spx_timestep_index_table){
  	delete [] this->spx_timestep_index_table;
  }
}

//----------------------------------------------------------------------------
int vtkMFIXReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro( << "Reading MFIX file");
  
  this->MakeMesh(output);

  return 1;
}

//----------------------------------------------------------------------------
void vtkMFIXReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << "\n";

  os << indent << "Number Of Nodes: " << this->NumberOfPoints << endl;
  os << indent << "Number Of Cells: " << this->NumberOfCells << endl;
  os << indent << "Number Of Cell Fields: " << this->NumberOfCellFields << endl;
}

//----------------------------------------------------------------------------
void vtkMFIXReader::MakeMesh(vtkUnstructuredGrid *output)
{
	output->Allocate();

	if(this->MakeMeshFlag == 0) {
		points->SetNumberOfPoints((imax2+1)*(jmax2+1)*(kmax2+1));
		
		//
		//  Cartesian type mesh
		//
		
		if( !strcmp(coordinates,"CARTESIAN") || (kmax2 == 1)) {
			double px = 0.0;
			double py = 0.0;
			double pz = 0.0;
			int cnt = 0;
			for (int k=0; k<= kmax2; k++){
				for (int j=0; j<= jmax2; j++){
					for (int i=0; i<= imax2; i++){
						points->InsertPoint(cnt, px, py, pz );
						cnt++;
						if( i == imax2 ) {
							px = px + DX->GetValue(i-1);
						} else { 
							px = px + DX->GetValue(i);
						}
					}
					px = 0.0;
					if ( j == jmax2) {
						py = py + DY->GetValue(j-1);
					} else { 
						py = py + DY->GetValue(j);
					}
				}
				py = 0.0;
				if ( k == kmax2) {
					pz = pz + DZ->GetValue(k-1);
				} else { 
					pz = pz + DZ->GetValue(k);
				}
			}
		} else {
		//
		//  Cylindrical Type Mesh
		//
			double px = 0.0;
			double py = 0.0;
			double pz = 0.0;
			double rx = 0.0;
			double ry = 0.0;
			double rz = 0.0;
			int cnt = 0;
			for (int k=0; k<= kmax2; k++){
				for (int j=0; j<= jmax2; j++){
					for (int i=0; i<= imax2; i++){
						points->InsertPoint(cnt, rx, ry, rz );
						cnt++;
						if( i == imax2 ) {
							px = px + DX->GetValue(i-1);
						} else if ( i == 0 ) {
							px = 0;
						} else { 
							px = px + DX->GetValue(i);
						}
						rx = px * cos(pz);
						rz = px * sin(pz) * -1;
					}
					px = 0.0;
					rx = 0.0;
					rz = 0.0;
					if ( j == jmax2) {
						py = py + DY->GetValue(j-1);
					} else { 
						py = py + DY->GetValue(j);
					}
					ry = py;
				}
				py = 0.0;
				ry = 0.0;
				if ( k == kmax2) {
					pz = pz + DZ->GetValue(k-1);
				} else { 
					pz = pz + DZ->GetValue(k);
				}
			}
		}
	

	//
	//  Let's put the points in a mesh
	//

		mesh->SetPoints(points);
	
		int p0 = 0;
		int count = 0;
		for (int k=0; k< kmax2; k++){
			for (int j=0; j< jmax2; j++){
				for (int i=0; i< imax2; i++){
					if ( FLAG->GetValue(count) < 10 ) {
						if ( !strcmp(coordinates,"CYLINDRICAL" ) && (kmax2 != 1) && ( i==0)) {
							aWedge->GetPointIds()->SetId( 0, p0-1);
							aWedge->GetPointIds()->SetId( 1, p0+1+((imax2+1)*(jmax2+1)));
							aWedge->GetPointIds()->SetId( 2, p0+1);
							aWedge->GetPointIds()->SetId( 3, p0+imax2);
							aWedge->GetPointIds()->SetId( 4, p0+2+imax2+((imax2+1)*(jmax2+1)));
							aWedge->GetPointIds()->SetId( 5, p0+2+imax2);
							mesh->InsertNextCell(aHexahedron->GetCellType(), aWedge->GetPointIds());
						} else {
							aHexahedron->GetPointIds()->SetId( 0, p0);
							aHexahedron->GetPointIds()->SetId( 1, p0+1);
							aHexahedron->GetPointIds()->SetId( 2, p0+1+((imax2+1)*(jmax2+1)));
							aHexahedron->GetPointIds()->SetId( 3, p0+((imax2+1)*(jmax2+1)));
							aHexahedron->GetPointIds()->SetId( 4, p0+1+imax2);
							aHexahedron->GetPointIds()->SetId( 5, p0+2+imax2);
							aHexahedron->GetPointIds()->SetId( 6, p0+2+imax2+((imax2+1)*(jmax2+1)));
							aHexahedron->GetPointIds()->SetId( 7, p0+1+imax2+((imax2+1)*(jmax2+1)));
							mesh->InsertNextCell(aHexahedron->GetCellType(), aHexahedron->GetPointIds());
						}
					}
					p0++;
					count++;
				}
				p0++;
			}
			p0 = p0 + imax2+1;
		}
		
 		cell_data_array = new vtkFloatArray * [variable_names->GetMaxId()+2];
		for (int j=0; j<=variable_names->GetMaxId(); j++){
			cell_data_array[ j ] = vtkFloatArray::New();
			cell_data_array[ j ]->SetName(variable_names->GetValue(j));
			cell_data_array[ j ]->SetNumberOfComponents(variable_components->GetValue(j));
		}
		
		min = new float [variable_names->GetMaxId()+1];
		max = new float [variable_names->GetMaxId()+1];
		veclen = new int [variable_names->GetMaxId()+1];

		this->MakeMeshFlag = 1;
	}
	
	output->DeepCopy(mesh);  // If mesh has already been made copy it to output

	int first = 0;
	for (int j=0; j<=variable_names->GetMaxId(); j++){
		
		if( this->CellDataArraySelection->GetArraySetting(j) == 1 ){
			if(variable_components->GetValue(j) == 1){
				GetVariableAtTimestep( j, this->TimeStep, cell_data_array[j]);
			} else {
				FillVectorVariable( j-3, j-2, j-1, cell_data_array[j]);
			}
			
			if(first == 0) {
				output->GetCellData()->SetScalars(cell_data_array[j]);
			} else {
				output->GetCellData()->AddArray(cell_data_array[j]);
			}
			
			double tempRange[2];
			cell_data_array[j]->GetRange(tempRange, j);
			min[j] = tempRange[0];
			max[j] = tempRange[1];
			veclen[j] = 1;
			first = 1;
		}
	}
}


//----------------------------------------------------------------------------
int vtkMFIXReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  
  if(this->RequestInformationFlag == 0) {
	
	if ( !this->FileName )
	{
	this->NumberOfPoints = 0;
	this->NumberOfCells = 0;
	
	vtkErrorMacro("No filename specified");
	return 0;
	}
	
	SetProjectName(this->FileName);
	ReadRes0();
	
	CreateVariableNames();
	GetTimeSteps();
	CalculateMaxTimeStep();
	MakeTimeStepTable(variable_names->GetMaxId()+1);
	GetNumberOfVariablesInSPXFiles();
	MakeSPXTimeStepIndexTable(variable_names->GetMaxId()+1);
	for (int j=0; j<=variable_names->GetMaxId(); j++){
		this->CellDataArraySelection->AddArray(variable_names->GetValue(j));
	}
	
	this->NumberOfPoints = (imax2+1)*(jmax2+1)*(kmax2+1);
	this->NumberOfCells = ijkmax2;
	this->NumberOfCellFields = variable_names->GetMaxId()+1;
	this->NumberOfTimeSteps = this->max_timestep;
	this->TimeStepRange[0] = 0;  
	this->TimeStepRange[1] = this->NumberOfTimeSteps-1;
	this->RequestInformationFlag = 1;
	
  }
  
  return 1;
}
//----------------------------------------------------------------------------
int vtkMFIXReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkMFIXReader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkMFIXReader::GetCellArrayStatus(const char* name)
{              
  return this->CellDataArraySelection->ArrayIsEnabled(name);

}
  
//----------------------------------------------------------------------------
void vtkMFIXReader::SetCellArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->CellDataArraySelection->EnableArray(name);
    }
  else
    {
    this->CellDataArraySelection->DisableArray(name);
    }

}

//----------------------------------------------------------------------------
void vtkMFIXReader::DisableAllCellArrays()
{
    this->CellDataArraySelection->DisableAllArrays();

}

//----------------------------------------------------------------------------
void vtkMFIXReader::EnableAllCellArrays()
{
    this->CellDataArraySelection->EnableAllArrays();

}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetCellDataRange(int cellComp, int index, float *min, float *max)
{
  if (index >= this->veclen[cellComp] || index < 0)
    {
    index = 0;  // if wrong index, set it to zero
    }
  *min = this->min[cellComp];
  *max = this->max[cellComp];

}

void vtkMFIXReader::SetProjectName (char *infile) {

	char temp[120];
	
	int len = strlen(infile);
	
	strncpy(run_name, infile, len-4);
	
}

void vtkMFIXReader::RestartVersionNumber(char* buffer)
{
   char s1[120];
   char s2[120];
   
   sscanf(buffer,"%s %s %f", s1, s2, &version_number);
   
   strncpy(version, buffer, 100);
   
}

void vtkMFIXReader::GetInt(istream& in, int &val)
{
	in.read( (char*)&val,sizeof(int));
	SWAP_INT(val);
}

void vtkMFIXReader::SWAP_INT(int &value)
{
	static char Swapped[4];
	
	int * Addr = &value;
	
	Swapped[0]=*((char*)Addr+3);
	Swapped[1]=*((char*)Addr+2);
	Swapped[2]=*((char*)Addr+1);
	Swapped[3]=*((char*)Addr  );
	
	value = *(reinterpret_cast<int*>(Swapped));
}

void vtkMFIXReader::SWAP_DOUBLE(double &value)
{
	static char Swapped[8];
	
	double * Addr = &value;
	
	Swapped[0]=*((char*)Addr+7);
	Swapped[1]=*((char*)Addr+6);
	Swapped[2]=*((char*)Addr+5);
	Swapped[3]=*((char*)Addr+4);
	Swapped[4]=*((char*)Addr+3);
	Swapped[5]=*((char*)Addr+2);
	Swapped[6]=*((char*)Addr+1);
	Swapped[7]=*((char*)Addr  );
	
	value = *(reinterpret_cast<double*>(Swapped));
}

void vtkMFIXReader::SWAP_FLOAT(float &value)
{
	static char Swapped[4];
	
	float * Addr = &value;
	
	Swapped[0]=*((char*)Addr+3);
	Swapped[1]=*((char*)Addr+2);
	Swapped[2]=*((char*)Addr+1);
	Swapped[3]=*((char*)Addr  );
	
	value = *(reinterpret_cast<float*>(Swapped));
}

void vtkMFIXReader::GetDouble(istream& in, double& val)
{
	in.read( (char*)&val,sizeof(double));
	SWAP_DOUBLE(val);
}

void vtkMFIXReader::SkipBytes(istream& in, int n)
{
	in.read(buffer,n); // maybe seekg instead
}

void vtkMFIXReader::IN_BIN_512(istream& in, vtkDoubleArray *v, int n)
{
   const int nr = 512/sizeof(double);

   double array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( (char*)&array , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
	      double temp = array[j];
              SWAP_DOUBLE(temp);
              v->InsertValue( c, temp);
              ++c;
           }
       }
   }
}

void vtkMFIXReader::IN_BIN_512I(istream& in, vtkIntArray *v, int n)
{
   const int nr = 512/sizeof(int);

   int array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( (char*)&array , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
	      int temp = array[j];
              SWAP_INT(temp);
              v->InsertValue( c, temp);
              ++c;
           }
       }
   }
}

void vtkMFIXReader::IN_BIN_512R(istream& in, vtkFloatArray *v, int n)
{
   const int nr = 512/sizeof(float);
   float array[nr];

   int num_records;
   
   if (!in) {
      cout << "Error opening file" << endl;
   }
  
   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;
   
   int c = 0;
   int cnt = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( (char*)&array , 512 );
       
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              float temp = array[j];
              SWAP_FLOAT(temp);
	      if( FLAG->GetValue(c) < 10) {
	        v->InsertValue(cnt, temp);
	        cnt++;
	      }
              ++c;
           }
       }
   }
}

void vtkMFIXReader::ReadRes0()
{
   int i,l,n,lc;
   
   int DIMENSION_USR = 5;

   ifstream in(this->FileName,ios::binary);

   if (!in)
   {
      cout << "could not open file" << endl;
      return;
   }

   buffer[512] = '\0';

   // version : record 1
   memset(buffer,0,513);
   in.read(buffer,512);

   RestartVersionNumber(buffer);
   

   // skip 2 linesline : records 2 and 3
   in.read(buffer,512);
   in.read(buffer,512);

   // imin1 etc : record 4
   memset(buffer,0,513);
   
   if (version == "RES = 01.00")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 15 ints ... 4 doubles = 92 bytes
      SkipBytes(in,420);
   }
   else if(version == "RES = 01.01" || version == "RES = 01.02")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC);    GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 4 doubles = 100 bytes
      SkipBytes(in,412);
   }
   else if(version == "RES = 01.03")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 5 doubles = 108 bytes
      SkipBytes(in,404);
   }
   else if(version == "RES = 01.04")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 18 ints ... 5 doubles = 112 bytes
      SkipBytes(in,400);
   }
   else if(version == "RES = 01.05")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 19 ints ... 5 doubles = 116 bytes
      SkipBytes(in,396);
   }
   else
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
      GetDouble(in,C_e); GetDouble(in,C_f); GetDouble(in,Phi); GetDouble(in,Phi_w);
      // 19 ints ... 9 doubles = 148 bytes
      SkipBytes(in,364);
   }

   const int nr = 512/sizeof(float);

   if ( ijkmax2%nr == 0)
      spx_records_per_timestep = ijkmax2/nr;
   else
      spx_records_per_timestep = 1 + ijkmax2/nr;

   
   // C , C_name and nmax

   NMAX->Resize(MMAX+1);
   for (int lc=0; lc<MMAX+1; ++lc){
	NMAX->InsertValue(lc, 1);
   }
   
   C->Resize(DIM_C);

   if (version_number > 1.04)
   {
      IN_BIN_512 (in, C, DIM_C);

      for (lc=0; lc<DIM_C; ++lc) in.read(buffer,512);  // c_name[] 
 
      if (version_number < 1.12)
          IN_BIN_512I(in, NMAX,MMAX+1);
      else
      {
          // what is the diff between this and above ??? 
          for (lc=0; lc<MMAX+1; ++lc) {
	      int temp;
              GetInt(in,temp);
	      NMAX->InsertValue(lc, temp);
	  }

          SkipBytes(in,512-(MMAX+1)*sizeof(int));
      }
   }
  
   DX->Resize(imax2);
   DY->Resize(jmax2);
   DZ->Resize(kmax2);

   IN_BIN_512(in, DX,imax2);
   IN_BIN_512(in, DY,jmax2);
   IN_BIN_512(in, DZ,kmax2);
	
   // run_name etc.
   
   memset(units,0,17);
   memset(coordinates,0,17);
   
   in.read(buffer,120);      // run_name , description
   in.read(units,16);        // units
   in.read(buffer,16);       // run_type
   in.read(coordinates,16);  // coordinates 
   
   SkipBytes(in,512-168);

   char tmp[17];
   
   memset(tmp,0,17);

   int ic = 0;
   for (i=0; i<17; ++i)
   {
        if (units[i] != ' ') tmp[ic++] = units[i];
   }

   memset(tmp,0,17);

   ic = 0;
   for (i=0; i<17; ++i)
   {
        if (coordinates[i] != ' ') tmp[ic++] = coordinates[i];
   }
   strcpy(coordinates,tmp);
   
   if (version_number >= 1.04)
   {
      tmpD->Resize(NMAX->GetValue(0));
      IN_BIN_512(in, tmpD, NMAX->GetValue(0));             // MW_g
      for (i=0; i<MMAX; ++i) in.read(buffer,512);  // MW_s
   }
   in.read(buffer,512);  // D_p etc.

   // read in the "DIM_IC" variables (and ignore ... not used by ani_mfix)
   tmpI->Resize(DIM_IC);
   tmpD->Resize(DIM_IC);

   IN_BIN_512(in, tmpD,DIM_IC);  // ic_x_w
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_x_e
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_y_s
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_y_n
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_z_b
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_z_t
   
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_i_w
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_i_e
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_j_s
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_j_n
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_k_b
   IN_BIN_512I(in, tmpI,DIM_IC);  // ic_k_t
   
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_ep_g
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_p_g
   IN_BIN_512(in, tmpD,DIM_IC);  // ic_t_g

   if (version_number < 1.15)
   {
      IN_BIN_512(in,tmpD,DIM_IC);  // ic_t_s(1,1)
      IN_BIN_512(in,tmpD,DIM_IC);  // ic_t_s(1,2) or ic_tmp 
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX->GetValue(0); ++i) IN_BIN_512(in,tmpD,DIM_IC); // ic_x_g
   }

   IN_BIN_512(in,tmpD,DIM_IC); // ic_u_g
   IN_BIN_512(in,tmpD,DIM_IC); // ic_v_g
   IN_BIN_512(in,tmpD,DIM_IC); // ic_w_g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,tmpD,DIM_IC); // ic_rop_s
      IN_BIN_512(in,tmpD,DIM_IC); // ic_u_s
      IN_BIN_512(in,tmpD,DIM_IC); // ic_v_s
      IN_BIN_512(in,tmpD,DIM_IC); // ic_w_s
      
      if (version_number >= 1.15)
      {
         IN_BIN_512(in,tmpD,DIM_IC); // ic_t_s
      }
      
      if (version_number >= 1.04)
      {
         for (n=0; n<NMAX->GetValue(lc+1); ++n)
            IN_BIN_512(in,tmpD,DIM_IC); // ic_x_s
      }
   }

   // read in the "DIM_BC" variables (and ignore ... not used by ani_mfix)
   tmpI->Resize(DIM_BC);
   tmpD->Resize(DIM_BC);

   IN_BIN_512(in,tmpD,DIM_BC); // bc_x_w
   IN_BIN_512(in,tmpD,DIM_BC); // bc_x_e
   IN_BIN_512(in,tmpD,DIM_BC); // bc y s
   IN_BIN_512(in,tmpD,DIM_BC); // bc y n
   IN_BIN_512(in,tmpD,DIM_BC); // bc z b
   IN_BIN_512(in,tmpD,DIM_BC);  // bc z t
   IN_BIN_512I(in,tmpI,DIM_BC);  // bc i w
   IN_BIN_512I(in,tmpI,DIM_BC); // bc i e
   IN_BIN_512I(in,tmpI,DIM_BC); // bc j s
   IN_BIN_512I(in,tmpI,DIM_BC); // bc j n
   IN_BIN_512I(in,tmpI,DIM_BC); // bc k b
   IN_BIN_512I(in,tmpI,DIM_BC); // bc k t
   IN_BIN_512(in,tmpD,DIM_BC); // bc ep g
   IN_BIN_512(in,tmpD,DIM_BC); // bc p g
   IN_BIN_512(in,tmpD,DIM_BC); // bc t g

   if (version_number < 1.15)
   {
      IN_BIN_512(in,tmpD,DIM_BC); // bc_t_s(1,1)
      IN_BIN_512(in,tmpD,DIM_BC); // bc_t_s(1,1) or bc_tmp
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX->GetValue(0); ++i) IN_BIN_512(in,tmpD,DIM_BC); // bc_x_g
   }

   IN_BIN_512(in,tmpD,DIM_BC); // bc u g
   IN_BIN_512(in,tmpD,DIM_BC); // bc v g
   IN_BIN_512(in,tmpD,DIM_BC); // bc w g
   IN_BIN_512(in,tmpD,DIM_BC); // bc ro g
   IN_BIN_512(in,tmpD,DIM_BC); // bc_rop_g
   IN_BIN_512(in,tmpD,DIM_BC); // bc volflow g
   IN_BIN_512(in,tmpD,DIM_BC); // bc massflow g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,tmpD,DIM_BC); // bc rop s
      IN_BIN_512(in,tmpD,DIM_BC); // bc u s
      IN_BIN_512(in,tmpD,DIM_BC); // bc v s
      
      if (version_number >= 1.04)
      {
         IN_BIN_512(in,tmpD,DIM_BC); // bc w s

         if (version_number >= 1.15)
         {
            IN_BIN_512(in,tmpD,DIM_BC); // bc t s
         }
         for (n=0; n<NMAX->GetValue(lc+1); ++n)
         {      
            IN_BIN_512(in,tmpD,DIM_BC); // bc x s
         }
      }
      IN_BIN_512(in,tmpD,DIM_BC); // bc volflow s
      IN_BIN_512(in,tmpD,DIM_BC); // bc massflow s
   }


   if (version == "RES = 01.00")
      l = 10;
   else
      l = DIM_BC;

   for (lc=0; lc<l; ++lc) in.read(buffer,512); // BC TYPE

   FLAG->Resize(ijkmax2);
   IN_BIN_512I(in, FLAG,ijkmax2);

   // DIM_IS varibles (not needed by ani_mfix)
   tmpI->Resize(DIM_IS);
   tmpD->Resize(DIM_IS);

   if (version_number >= 1.04)
   {
      IN_BIN_512(in,tmpD,DIM_IS); // is x w
      IN_BIN_512(in,tmpD,DIM_IS); // is x e
      IN_BIN_512(in,tmpD,DIM_IS); // is y s
      IN_BIN_512(in,tmpD,DIM_IS); // is y n
      IN_BIN_512(in,tmpD,DIM_IS); // is z b
      IN_BIN_512(in,tmpD,DIM_IS); // is z t
      IN_BIN_512I(in,tmpI,DIM_IS); // is i w
      IN_BIN_512I(in,tmpI,DIM_IS); // is i e
      IN_BIN_512I(in,tmpI,DIM_IS); // is j s
      IN_BIN_512I(in,tmpI,DIM_IS); // is j n
      IN_BIN_512I(in,tmpI,DIM_IS); // is k b
      IN_BIN_512I(in,tmpI,DIM_IS); // is k t
      IN_BIN_512(in,tmpD,DIM_IS);  // is_pc(1,1)
      IN_BIN_512(in,tmpD,DIM_IS);  // is_pc(1,2)
     
      if (version_number >= 1.07)
      {
         for (l=0; l<MMAX; ++l) IN_BIN_512(in,tmpD,DIM_IS); // is_vel_s
      }

      for (lc=0; lc<DIM_IS; ++lc) in.read(buffer,512); // is_type
   }

   if (version_number >= 1.08) in.read(buffer,512);
   
   
   if (version_number >= 1.09) 
   {
      in.read(buffer,512);
      
      if (version_number >= 1.5)
      {
         GetInt(in,nspx_use);
         SkipBytes(in,508);
      }
      
      for (lc=0; lc< nspx_use; ++lc) in.read(buffer,512); // spx_dt
      
      for (lc=0; lc<MMAX+1; ++lc) in.read(buffer,512);    // species_eq
      
      tmpD->Resize(DIMENSION_USR);
      
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr_dt
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr x w
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr x e
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr y s
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr y n
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr z b
      IN_BIN_512(in,tmpD,DIMENSION_USR); // usr z t
     
      for (lc=0; lc<DIMENSION_USR; ++lc) in.read(buffer,512);    // usr_ext etc.
      
          
      tmpD->Resize(DIM_IC);      
      IN_BIN_512(in,tmpD,DIM_IC); // ic_p_star
      IN_BIN_512(in,tmpD,DIM_IC); // ic_l_scale
      for (lc=0; lc<DIM_IC; ++lc) in.read(buffer,512);    // ic_type
          
      tmpD->Resize(DIM_BC);      
      IN_BIN_512(in,tmpD,DIM_BC); // bc_dt_0
      IN_BIN_512(in,tmpD,DIM_BC); // bc_jet_g0
      IN_BIN_512(in,tmpD,DIM_BC); // bc_dt_h
      IN_BIN_512(in,tmpD,DIM_BC); // bc_jet_gh
      IN_BIN_512(in,tmpD,DIM_BC); // bc_dt_l
      IN_BIN_512(in,tmpD,DIM_BC); // bc_jet_gl
      }
      
      
      if (version_number >= 1.1)  in.read(buffer,512);  // mu_gmax
      if (version_number >= 1.11) in.read(buffer,512);  // x_ex , model_b
      
      if (version_number >= 1.12)
      {
         in.read(buffer,512);   // p_ref , etc.
         in.read(buffer,512);   // leq_it , leq_method
    
         IN_BIN_512(in,tmpD,DIM_BC); // bc_hw_g
         IN_BIN_512(in,tmpD,DIM_BC); // bc_uw_g
         IN_BIN_512(in,tmpD,DIM_BC); // bc_vw_g
         IN_BIN_512(in,tmpD,DIM_BC); // bc_ww_g
    
         for (lc=0; lc<MMAX; ++lc)
         {
            IN_BIN_512(in,tmpD,DIM_BC); // bc_hw_s
            IN_BIN_512(in,tmpD,DIM_BC); // bc_uw_s
            IN_BIN_512(in,tmpD,DIM_BC); // bc_vw_s
            IN_BIN_512(in,tmpD,DIM_BC); // bc_ww_s
         }
      }
      
      if (version_number >= 1.13) in.read(buffer,512);    // momentum_x_eq , etc.
      if (version_number >= 1.14) in.read(buffer,512);    // detect_small
      
      if (version_number >= 1.15)
      {
         in.read(buffer,512);    // k_g0 , etc.
     
         tmpD->Resize(DIM_IC);   
     
         IN_BIN_512(in,tmpD,DIM_IC); // ic_gama_rg
         IN_BIN_512(in,tmpD,DIM_IC); // ic_t_rg
        
         for (lc=0; lc<MMAX; ++lc)
         {
            IN_BIN_512(in,tmpD,DIM_IC); // ic_gama_rs
            IN_BIN_512(in,tmpD,DIM_IC); // ic_t_rs
         }
     }
     
     if (version_number >= 1.2) in.read(buffer,512); // norm_g , norm_s
 
     if (version_number >= 1.3)
     {
        GetInt(in,NScalar);
        SkipBytes(in,sizeof(double)); // tol_resid_scalar

        int DIM_tmp;
        GetInt(in,DIM_tmp);
        SkipBytes(in,512-sizeof(double)-2*sizeof(int));
    
        tmpI->Resize(DIM_tmp);
        IN_BIN_512I(in,tmpI,DIM_tmp);  // Phase4Scalar;
     }
        
     if (version_number >= 1.5)
     {
        GetInt(in,nRR);
        SkipBytes(in,508);
     }
     
     if (version_number >= 1.5999)
     {
        int tmp;
        GetInt(in,tmp);
        SkipBytes(in,508);
    
        if (tmp != 0) bKepsilon = true;
     }      
   

}

void vtkMFIXReader::CreateVariableNames()
{

    char fname[256];
    
    int cnt = 0;
    
    for (int i=0; i<nspx_use; ++i)
    {
    
	for(int k=0;k<sizeof(fname);k++){
		fname[k]=0;
	}
	strncpy(fname, this->FileName, strlen(this->FileName)-4);
	
	if(i==0){
		strcat(fname, ".SP1");
	} else if (i==1) {
		strcat(fname, ".SP2");
	} else if (i==2) {
		strcat(fname, ".SP3");
	} else if (i==3) {
		strcat(fname, ".SP4");
	} else if (i==4) {
		strcat(fname, ".SP5");
	} else if (i==5) {
		strcat(fname, ".SP6");
	} else if (i==6) {
		strcat(fname, ".SP7");
	} else if (i==7) {
		strcat(fname, ".SP8");
	} else if (i==8) {
		strcat(fname, ".SP9");
	} else if (i==9) {
		strcat(fname, ".SPA");
	} else {
		strcat(fname, ".SPB");
	} 
	
	ifstream in(fname,ios::binary);
        if (in) // file exists
        {
	    this->SpxFileExists->InsertValue(i, 1);
	
            switch (i+1)
            {

            case 1:
                {
                    variable_names->InsertValue(cnt++,"EP_g");
		    variableIndexToSPX->InsertValue(cnt-1, 1);
		    variable_components->InsertValue(cnt-1, 1);
                    break;
                }

            case 2:
                {
                    variable_names->InsertValue(cnt++,"P_g");
		    variableIndexToSPX->InsertValue(cnt-1, 2);
		    variable_components->InsertValue(cnt-1, 1);
                    variable_names->InsertValue(cnt++,"P_star");
		    variableIndexToSPX->InsertValue(cnt-1, 2);
		    variable_components->InsertValue(cnt-1, 1);
                    break;
                }

            case 3:
                {
                    variable_names->InsertValue(cnt++,"U_g");
		    variableIndexToSPX->InsertValue(cnt-1, 3);
		    variable_components->InsertValue(cnt-1, 1);
                    variable_names->InsertValue(cnt++,"V_g");
		    variableIndexToSPX->InsertValue(cnt-1, 3);
		    variable_components->InsertValue(cnt-1, 1);
                    variable_names->InsertValue(cnt++,"W_g");
		    variableIndexToSPX->InsertValue(cnt-1, 3);
		    variable_components->InsertValue(cnt-1, 1);
                    variable_names->InsertValue(cnt++,"Gas Velocity");
		    variableIndexToSPX->InsertValue(cnt-1, 3);
		    variable_components->InsertValue(cnt-1, 3);
                    break;
                }


            case 4:
                {
		    char us[120];
		    char vs[120];
		    char ws[120];
		    char sv[120];
		    char temp[120];
		    
		    
                    for (int i=0; i<MMAX; ++i)
                    {
		    	for(int k=0;k<sizeof(us);k++) {
			   us[k]=0;
			}	
		    	for(int k=0;k<sizeof(vs);k++) {
			   vs[k]=0;
			}	
		    	for(int k=0;k<sizeof(ws);k++) {
			   ws[k]=0;
			}	
		    	for(int k=0;k<sizeof(sv);k++) {
			   sv[k]=0;
			}	
		        strcpy(us, "U_s_");
		        strcpy(vs, "V_s_");
		        strcpy(ws, "W_s_");
		        strcpy(sv, "Solids_Velocity_");
		    	sprintf(temp, "%d", i+1);
			strcat(us, temp);
			strcat(vs, temp);
			strcat(ws, temp);
			strcat(sv, temp);
                        variable_names->InsertValue(cnt++, us);
		    	variableIndexToSPX->InsertValue(cnt-1, 4);
			variable_components->InsertValue(cnt-1, 1);
                        
			variable_names->InsertValue(cnt++, vs);
		    	variableIndexToSPX->InsertValue(cnt-1, 4);
			variable_components->InsertValue(cnt-1, 1);
                        
			variable_names->InsertValue(cnt++, ws);
		    	variableIndexToSPX->InsertValue(cnt-1, 4);
			variable_components->InsertValue(cnt-1, 1);
                        
			variable_names->InsertValue(cnt++, sv);
		    	variableIndexToSPX->InsertValue(cnt-1, 4);
			variable_components->InsertValue(cnt-1, 3);

                    }

                    break;
                }

            case 5:
                {
		    char rops[120];
		    char temp[120];

                    for (int i=0; i<MMAX; ++i)
                    {
		    	for(int k=0;k<sizeof(rops);k++) {
			   rops[k]=0;
			}	
			strcpy(rops, "ROP_s_");
		    	sprintf(temp, "%d", i+1);
			strcat(rops, temp);
                        variable_names->InsertValue(cnt++, rops);
		        variableIndexToSPX->InsertValue(cnt-1, 5);
			variable_components->InsertValue(cnt-1, 1);
                    }

                    break;
                }

            case 6:
                {
                    variable_names->InsertValue(cnt++, "T_g");
	            variableIndexToSPX->InsertValue(cnt-1, 6);
		    variable_components->InsertValue(cnt-1, 1);

                    if (version_number <= 1.15)
                    {
                        variable_names->InsertValue(cnt++, "T_s_1");
		        variableIndexToSPX->InsertValue(cnt-1, 6);
			variable_components->InsertValue(cnt-1, 1);

                        if (MMAX > 1){ 
                            variable_names->InsertValue(cnt++, "T_s_2");
			    variableIndexToSPX->InsertValue(cnt-1, 6);
			    variable_components->InsertValue(cnt-1, 1);
			} else {
                            variable_names->InsertValue(cnt++, "T_s_2_not_used");
			    variableIndexToSPX->InsertValue(cnt-1, 6);
			    variable_components->InsertValue(cnt-1, 1);
			}
                    }
                    else
                    {
			char ts[120];
			char temp[120];
			
			
			for (int i=0; i<MMAX; ++i)
			{
				for(int k=0;k<sizeof(ts);k++) {
				ts[k]=0;
				}	
				strcpy(ts, "T_s_");
				sprintf(temp, "%d", i+1);
				strcat(ts, temp);
				variable_names->InsertValue(cnt++, ts);
			        variableIndexToSPX->InsertValue(cnt-1, 6);
			        variable_components->InsertValue(cnt-1, 1);
			}
                    }

                    break;
                }

            case 7:
                {
			char var[120];
			char temp[120];

                    for (int i=0; i<NMAX->GetValue(0); ++i)
                    {
			for(int k=0;k<sizeof(var);k++) {
			   var[k]=0;
			}	
			strcpy(var, "X_g_");
			sprintf(temp, "%d", i+1);
			strcat(var, temp);
			variable_names->InsertValue(cnt++, var);
			variableIndexToSPX->InsertValue(cnt-1, 7);
			variable_components->InsertValue(cnt-1, 1);
                    }


                    for (int m=0; m<MMAX; ++m)
                    {
                        for (int i=0; i<NMAX->GetValue(m); ++i)
                        {
			    char temp1[120];
			    char temp2[120];
			    for(int k=0;k<sizeof(var);k++) {
				var[k]=0;
			    }	
			    strcpy(var, "X_s_");
			    sprintf(temp1, "%d", m+1);
			    sprintf(temp2, "%d", i+1);
			    strcat(var, temp1);
			    strcat(var, "_");
			    strcat(var, temp2);
			    variable_names->InsertValue(cnt++, var);
    			    variableIndexToSPX->InsertValue(cnt-1, 7);
			    variable_components->InsertValue(cnt-1, 1);
                        }
                    }


                    break;
                }


            case 8:
                {
		    char var[120];
		    char temp[120];
                    for (int i=0; i<MMAX; ++i)
                    {
			for(int k=0;k<sizeof(var);k++) {
			  var[k]=0;
			}	
			strcpy(var, "Theta_m_");
			sprintf(temp, "%d", i+1);
			strcat(var, temp);
			variable_names->InsertValue(cnt++, var);
  		        variableIndexToSPX->InsertValue(cnt-1, 8);
			variable_components->InsertValue(cnt-1, 1);
                    }

                    break;
                }


            case 9:
                {
		    char var[120];
		    char temp[120];

                    for (int i=0; i<NScalar; ++i)
                    {
			for(int k=0;k<sizeof(var);k++) {
			   var[k]=0;
			}	
		        strcpy(var, "Scalar_");
			sprintf(temp, "%d", i+1);
			strcat(var, temp);
			variable_names->InsertValue(cnt++, var);
			variableIndexToSPX->InsertValue(cnt-1, 9);
			variable_components->InsertValue(cnt-1, 1);
                    }

                    break;
                }


            case 10:
                {
		    char var[120];
		    char temp[120];

                    for (int i=0; i<nRR; ++i)
                    {
			for(int k=0;k<sizeof(var);k++) {
			   var[k]=0;
			}	
		        strcpy(var, "RRates_");
			sprintf(temp, "%d", i+1);
			strcat(var, temp);
			variable_names->InsertValue(cnt++, var);
		        variableIndexToSPX->InsertValue(cnt-1, 10);
			variable_components->InsertValue(cnt-1, 1);
                    }

                    break;
                }

            case 11:
                {
                    if (bKepsilon)
                    {
                        variable_names->InsertValue(cnt++, "k_turb_g");
			variableIndexToSPX->InsertValue(cnt-1, 11);
			variable_components->InsertValue(cnt-1, 1);
                        variable_names->InsertValue(cnt++, "e_turb_g");
			variableIndexToSPX->InsertValue(cnt-1, 11);
			variable_components->InsertValue(cnt-1, 1);
                    }

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }


            }


        } else {
	        this->SpxFileExists->InsertValue(i, 0);
	}
	
    }

}

void vtkMFIXReader::GetTimeSteps()
{
	int next_rec, num_rec;
	
	char fname[256];
 	int cnt = 0;
	
	for (int i=0; i<nspx_use; ++i)
	{
		
		for(int k=0;k<sizeof(fname);k++){
			fname[k]=0;
		}
		strncpy(fname, this->FileName, strlen(this->FileName)-4);
		
		if(i==0){
			strcat(fname, ".SP1");
		} else if (i==1) {
			strcat(fname, ".SP2");
		} else if (i==2) {
			strcat(fname, ".SP3");
		} else if (i==3) {
			strcat(fname, ".SP4");
		} else if (i==4) {
			strcat(fname, ".SP5");
		} else if (i==5) {
			strcat(fname, ".SP6");
		} else if (i==6) {
			strcat(fname, ".SP7");
		} else if (i==7) {
			strcat(fname, ".SP8");
		} else if (i==8) {
			strcat(fname, ".SP9");
		} else if (i==9) {
			strcat(fname, ".SPA");
		} else {
			strcat(fname, ".SPB");
		} 

		ifstream in(fname , ios::binary);

		int nvars=0;
		if (in) // file exists
		{
			in.clear();
			in.seekg( 1024, ios::beg ); 
			
			in.read( (char*)&next_rec,sizeof(int) );
			SWAP_INT(next_rec);
			
			in.read( (char*)&num_rec,sizeof(int) );
			SWAP_INT(num_rec);
			

			switch (i+1)
			{
				case 1: nvars = 1; break;
				case 2: nvars = 2; break; 
				case 3: nvars = 4; break;
				case 4: nvars = 4*MMAX; break;
				case 5: nvars = MMAX; break;
				case 6:
				{
					if (version_number <= 1.15)
						nvars = 3;
					else
						nvars = MMAX + 1;
					break;
				}
				case 7:
				{
					nvars = NMAX->GetValue(0);
					for (int m=0; m<MMAX; ++m) {
						nvars += NMAX->GetValue(m);
					}
					break;
				}
				case 8: nvars = MMAX; break;
				case 9: nvars = NScalar; break;
				case 10: nvars = nRR; break; 
				case 11:
				{
					if (bKepsilon) nvars = 2;
					break;
				}
			}
			
			for(int j=0; j<nvars; j++){
				var_timesteps->InsertValue(cnt, (next_rec-4)/num_rec);
				cnt++;
			}
		}
		
	}	
}

void vtkMFIXReader::MakeTimeStepTable(int nvars)
{
	var_timestep_table->SetNumberOfComponents(nvars);
	
	for(int i=0; i<nvars; i++){
		int ts_increment = max_timestep/var_timesteps->GetValue(i);
		int ts = 1;
		for (int j=0; j<max_timestep; j++){
			var_timestep_table->InsertComponent(j, i, ts);
			ts_increment--;
			if(ts_increment <= 0) {
				ts_increment = max_timestep/var_timesteps->GetValue(i);
				ts++;
			}
			if(ts > var_timesteps->GetValue(i)) {
				ts = var_timesteps->GetValue(i);
			}
		}
	}
}

void vtkMFIXReader::GetVariableAtTimestep(int vari , int tstep, vtkFloatArray *v)
{
    // This routine opens and closes the file for each request.
    // Maybe keep all SPX files open, and just perform relative
    // moves to get to the correct location in the file

    // get filename that vaiable # vari is located in

    // assumptions : there are <10 solid phases,
    // <10 scalars and <10 ReactionRates (need to change this)

    char vname[256];
    strcpy(vname, variable_names->GetValue(vari));
    
    int spx = variableIndexToSPX->GetValue(vari);
    
    char fname[256];
    for(int k=0;k<sizeof(fname);k++){
       fname[k]=0;
    }

    strncpy(fname, this->FileName, strlen(this->FileName)-4);

    if(spx==1){
	strcat(fname, ".SP1");
    } else if (spx==2) {
	strcat(fname, ".SP2");
    } else if (spx==3) {
	strcat(fname, ".SP3");
    } else if (spx==4) {
	strcat(fname, ".SP4");
    } else if (spx==5) {
	strcat(fname, ".SP5");
    } else if (spx==6) {
	strcat(fname, ".SP6");
    } else if (spx==7) {
	strcat(fname, ".SP7");
    } else if (spx==8) {
	strcat(fname, ".SP8");
    } else if (spx==9) {
	strcat(fname, ".SP9");
    } else if (spx==10) {
	strcat(fname, ".SPA");
    } else {
	strcat(fname, ".SPB");
    } 

    int ind = (vari*max_timestep) + tstep;
    
    int nBytesSkip = spx_timestep_index_table[ind];

    ifstream in(fname,ios::binary);

    in.seekg(nBytesSkip,ios::beg);

    IN_BIN_512R (in, v, ijkmax2);
    

}

void vtkMFIXReader::MakeSPXTimeStepIndexTable(int nvars)
{    
	int spx_timestep_index_table_size = nvars * max_timestep;
	spx_timestep_index_table = new int [spx_timestep_index_table_size];
	
	int timestep;
	int spx;
	int nvars_in_spx;
	for(int i=0; i<nvars; i++){
		for (int j=0; j<max_timestep; j++){
			timestep = (int) var_timestep_table->GetComponent(j, i);
			spx = variableIndexToSPX->GetValue(i);
			nvars_in_spx = spx_to_nvar_table->GetValue(spx);
			int skip = var_to_skip_table->GetValue(i);
			
			int index = (3*512) + (timestep-1)*((nvars_in_spx*spx_records_per_timestep*512)+512) + 512 + (skip*spx_records_per_timestep*512);
			
			int ind = (i*max_timestep) + j;
			spx_timestep_index_table[ind] = index;

		}
	}

}

void vtkMFIXReader::CalculateMaxTimeStep()
{

//
//   What is the maximum timestep value?
//
	
  this->max_timestep = 0;
  for ( int i=0; i <= variable_names->GetMaxId() ; i++ ) {
	if(var_timesteps->GetValue(i) > this->max_timestep){
		this->max_timestep = var_timesteps->GetValue(i);
	}
  }

}

void vtkMFIXReader::GetNumberOfVariablesInSPXFiles()
{
  
  //
  //  How many variables are in each spx?
  //
  int spx_nvars = 0;
  int skip = 0;
  for(int j=1; j<nspx_use; j++) {
	for(int i=0;i<variable_names->GetMaxId()+1;i++){
		if((variableIndexToSPX->GetValue(i) == j) && (variable_components->GetValue(i) == 1)){
			spx_nvars++;
			var_to_skip_table->InsertValue(i,skip);
			skip++;
		}
	}
	spx_to_nvar_table->InsertValue(j, spx_nvars);
	spx_nvars = 0;
	skip = 0;
  }	
}

void vtkMFIXReader::FillVectorVariable( int xindex, int yindex, int zindex, vtkFloatArray *v)
{
	int range = cell_data_array[xindex]->GetMaxId();
	
	for(int i=0;i<=range;i++){
		v->InsertComponent(i, 0, cell_data_array[xindex]->GetValue(i));
		v->InsertComponent(i, 1, cell_data_array[yindex]->GetValue(i));
		v->InsertComponent(i, 2, cell_data_array[zindex]->GetValue(i));
	}
}
