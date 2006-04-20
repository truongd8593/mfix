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
// Please address all comments to Brian Dotson (brian.dotson@netl.doe.gov)
//

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
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

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
  this->Minimum = NULL;
  this->Maximum = NULL;
  this->VectorLength = NULL;
  this->CellDataArray = NULL;
  this->SPXTimestepIndexTable = NULL;
  this->DimensionIc = 5;
  this->DimensionBc = 5;
  this->DimensionC = 5;
  this->DimensionIs = 5;
  this->NumberOfSPXFilesUsed = 9;
  this->NumberOfScalars = 0;
  this->BkEpsilon = false;
  this->NumberOfReactionRates = 0;
  this->FileExtension[0] = '1';
  this->FileExtension[1] = '2';
  this->FileExtension[2] = '3';
  this->FileExtension[3] = '4';
  this->FileExtension[4] = '5';
  this->FileExtension[5] = '6';
  this->FileExtension[6] = '7';
  this->FileExtension[7] = '8';
  this->FileExtension[8] = '9';
  this->FileExtension[9] = 'A';
  this->FileExtension[10] = 'B';
  this->VersionNumber = 0;

  this->CellDataArray = NULL;
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->Points = vtkPoints::New();
  this->Mesh = vtkUnstructuredGrid::New();
  this->AHexahedron = vtkHexahedron::New();
  this->AWedge = vtkWedge::New();
  this->NMax = vtkIntArray::New();
  this->C = vtkDoubleArray::New();
  this->Dx = vtkDoubleArray::New();
  this->Dy = vtkDoubleArray::New();
  this->Dz = vtkDoubleArray::New();
  this->TempI = vtkIntArray::New();
  this->TempD = vtkDoubleArray::New();   
  this->Flag = vtkIntArray::New();
  this->VariableNames = vtkStringArray::New();
  this->VariableComponents = vtkIntArray::New();
  this->VariableIndexToSPX = vtkIntArray::New();
  this->VariableTimesteps = vtkIntArray::New();
  this->VariableTimestepTable = vtkIntArray::New();
  this->SPXToNVarTable = vtkIntArray::New();
  this->VariableToSkipTable = vtkIntArray::New();
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
  if ( this->FileName)
    {
    delete [] this->FileName;
    }

  for (int j = 0; j <= this->VariableNames->GetMaxId(); j++)
    {
    this->CellDataArray[j]->Delete();
    }

  this->CellDataArraySelection->Delete();
  this->Points->Delete();
  this->Mesh->Delete();
  this->AHexahedron->Delete();
  this->AWedge->Delete();
  this->NMax->Delete();
  this->C->Delete();
  this->Dx->Delete();
  this->Dy->Delete();
  this->Dz->Delete();
  this->TempI->Delete();
  this->TempD->Delete();
  this->Flag->Delete();
  this->VariableNames->Delete();
  this->VariableComponents->Delete();
  this->VariableIndexToSPX->Delete();
  this->VariableTimesteps->Delete();
  this->VariableTimestepTable->Delete();
  this->SPXToNVarTable->Delete();
  this->VariableToSkipTable->Delete();
  this->SpxFileExists->Delete();

  if (this->CellDataArray)
    {
    delete [] this->CellDataArray;
    }

  if (this->Minimum)
    {
    delete [] this->Minimum;
    }

  if (this->Maximum)
    {
    delete [] this->Maximum;
    }

  if (this->VectorLength)
    {
    delete [] this->VectorLength;
    }

  if (this->SPXTimestepIndexTable)
    {
    delete [] this->SPXTimestepIndexTable;
    }
}

//----------------------------------------------------------------------------
int vtkMFIXReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
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

  if (this->MakeMeshFlag == 0) 
    {
    Points->SetNumberOfPoints((IMaximum2+1)*(JMaximum2+1)*(KMaximum2+1));

    //
    //  Cartesian type mesh
    //
    if ( !strcmp(CoordinateSystem,"CARTESIAN") || (KMaximum2 == 1))
      {
      double px = 0.0;
      double py = 0.0;
      double pz = 0.0;
      int cnt = 0;
      for (int k = 0; k <= KMaximum2; k++)
        {
        for (int j = 0; j <= JMaximum2; j++)
          {
          for (int i = 0; i <= IMaximum2; i++)
            {
            Points->InsertPoint(cnt, px, py, pz );
            cnt++;
            if ( i == IMaximum2 )
              {
              px = px + Dx->GetValue(i-1);
              }
            else
              {
              px = px + Dx->GetValue(i);
              }
            }
          px = 0.0;
          if ( j == JMaximum2)
            {
            py = py + Dy->GetValue(j-1);
            }
          else
            {
            py = py + Dy->GetValue(j);
            }
          }
        py = 0.0;
        if ( k == KMaximum2)
          {
          pz = pz + Dz->GetValue(k-1);
          }
        else
          {
          pz = pz + Dz->GetValue(k);
          }
        }
      }
    else
      {
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
      for (int k = 0; k <= KMaximum2; k++)
        {
        for (int j = 0; j <= JMaximum2; j++)
          {
          for (int i = 0; i <= IMaximum2; i++)
            {
            Points->InsertPoint(cnt, rx, ry, rz );
            cnt++;
            if ( i == IMaximum2 )
              {
              px = px + Dx->GetValue(i-1);
              }
            else if ( i == 0 )
              {
              px = 0;
              }
            else
              {
              px = px + Dx->GetValue(i);
              }
            rx = px * cos(pz);
            rz = px * sin(pz) * -1;
            }
          px = 0.0;
          rx = 0.0;
          rz = 0.0;
          if ( j == JMaximum2)
            {
            py = py + Dy->GetValue(j-1);
            }
          else
            {
            py = py + Dy->GetValue(j);
            }
          ry = py;
          }
        py = 0.0;
        ry = 0.0;
        if ( k == KMaximum2)
          {
          pz = pz + Dz->GetValue(k-1);
          }
        else
          {
          pz = pz + Dz->GetValue(k);
          }
        }
      }

    //
    //  Let's put the points in a mesh
    //
    Mesh->SetPoints(Points);
    int p0 = 0;
    int count = 0;
    for (int k = 0; k < KMaximum2; k++)
      {
      for (int j = 0; j < JMaximum2; j++)
        {
        for (int i = 0; i < IMaximum2; i++)
          {
          if ( Flag->GetValue(count) < 10 )
            {
            if ( !strcmp(CoordinateSystem,"CYLINDRICAL" ) )
              {
              if (( k == (KMaximum2-2)) && (i != 1))
                {
                AHexahedron->GetPointIds()->SetId( 0, p0);
                AHexahedron->GetPointIds()->SetId( 1, p0+1);
                AHexahedron->GetPointIds()->SetId( 2,
                  (p0+1+((IMaximum2+1)*(JMaximum2+1)))-
                  ((IMaximum2+1)*(JMaximum2+1)*(KMaximum2-2)));
                AHexahedron->GetPointIds()->SetId( 3,
                  (p0+((IMaximum2+1)*(JMaximum2+1)))-
                  ((IMaximum2+1)*(JMaximum2+1)*(KMaximum2-2)));
                AHexahedron->GetPointIds()->SetId( 4, p0+1+IMaximum2);
                AHexahedron->GetPointIds()->SetId( 5, p0+2+IMaximum2);
                AHexahedron->GetPointIds()->SetId( 6, (p0+2+IMaximum2 +
                  ((IMaximum2+1)*(JMaximum2+1))) -
                  ((IMaximum2+1)*(JMaximum2+1) * (KMaximum2-2)));
                AHexahedron->GetPointIds()->SetId( 7, (p0+1+IMaximum2 +
                  ((IMaximum2+1)*(JMaximum2+1)))- 
                  ((IMaximum2+1)*(JMaximum2+1)*(KMaximum2-2)));
                Mesh->InsertNextCell(AHexahedron->GetCellType(), 
                  AHexahedron->GetPointIds());
                }
              else if ((k != (KMaximum2-2)) && (i != 1))
                {
                AHexahedron->GetPointIds()->SetId( 0, p0);
                AHexahedron->GetPointIds()->SetId( 1, p0+1);
                AHexahedron->GetPointIds()->SetId( 2, 
                  p0+1+((IMaximum2+1)*(JMaximum2+1)));
                AHexahedron->GetPointIds()->SetId( 3, 
                  p0+((IMaximum2+1)*(JMaximum2+1)));
                AHexahedron->GetPointIds()->SetId( 4, p0+1+IMaximum2);
                AHexahedron->GetPointIds()->SetId( 5, p0+2+IMaximum2);
                AHexahedron->GetPointIds()->SetId( 6, p0+2+IMaximum2+
                  ((IMaximum2+1)*(JMaximum2+1)));
                AHexahedron->GetPointIds()->SetId( 7, p0+1+IMaximum2+
                  ((IMaximum2+1)*(JMaximum2+1)));
                Mesh->InsertNextCell(AHexahedron->GetCellType(), 
                  AHexahedron->GetPointIds());
                }
              else if ( (k != (KMaximum2-2)) && (i == 1))
                {
                AWedge->GetPointIds()->SetId( 0, j*(IMaximum2+1));
                AWedge->GetPointIds()->SetId( 1, p0+1);
                AWedge->GetPointIds()->SetId( 2, p0+1+((IMaximum2+1)*
                  (JMaximum2+1)));
                AWedge->GetPointIds()->SetId( 3, (j+1)*(IMaximum2+1));
                AWedge->GetPointIds()->SetId( 4, p0+2+IMaximum2);
                AWedge->GetPointIds()->SetId( 5, p0+2+IMaximum2+
                  ((IMaximum2+1)*(JMaximum2+1)));
                Mesh->InsertNextCell(AWedge->GetCellType(), 
                  AWedge->GetPointIds());
                }
              else if (( k == (KMaximum2-2)) && (i == 1))
                {
                AWedge->GetPointIds()->SetId( 0, j*(IMaximum2+1));
                AWedge->GetPointIds()->SetId( 1, p0+1);
                AWedge->GetPointIds()->SetId( 2,
                 (p0+1+((IMaximum2+1)*(JMaximum2+1)))-((IMaximum2+1)
                 *(JMaximum2+1)*(KMaximum2-2)));
                AWedge->GetPointIds()->SetId( 3, (j+1)*(IMaximum2+1));
                AWedge->GetPointIds()->SetId( 4, p0+2+IMaximum2);
                AWedge->GetPointIds()->SetId( 5, (p0+2+IMaximum2 +
                  ((IMaximum2+1)*(JMaximum2+1))) -((IMaximum2+1)
                  *(JMaximum2+1)*(KMaximum2-2)));
                Mesh->InsertNextCell(AWedge->GetCellType(), 
                  AWedge->GetPointIds());
                }
              }
            else
              {
              AHexahedron->GetPointIds()->SetId( 0, p0);
              AHexahedron->GetPointIds()->SetId( 1, p0+1);
              AHexahedron->GetPointIds()->SetId( 2, p0+1+((IMaximum2+1)
                *(JMaximum2+1)));
              AHexahedron->GetPointIds()->SetId( 3, p0+((IMaximum2+1)
                *(JMaximum2+1)));
              AHexahedron->GetPointIds()->SetId( 4, p0+1+IMaximum2);
              AHexahedron->GetPointIds()->SetId( 5, p0+2+IMaximum2);
              AHexahedron->GetPointIds()->SetId( 6, p0+2+IMaximum2 +
                ((IMaximum2+1)*(JMaximum2+1)));
              AHexahedron->GetPointIds()->SetId( 7, p0+1+IMaximum2 + 
                ((IMaximum2+1)*(JMaximum2+1)));
              Mesh->InsertNextCell(AHexahedron->GetCellType(), 
                AHexahedron->GetPointIds());
              }
            }
          p0++;
          count++;
          }
        p0++;
        }
      p0 = p0 + IMaximum2+1;
      }

    CellDataArray = new vtkFloatArray * [VariableNames->GetMaxId()+2];
    for (int j = 0; j <= VariableNames->GetMaxId(); j++)
      {
      CellDataArray[ j ] = vtkFloatArray::New();
      CellDataArray[ j ]->SetName(VariableNames->GetValue(j));
      CellDataArray[ j ]->
        SetNumberOfComponents(VariableComponents->GetValue(j));
      }

    Minimum = new float [VariableNames->GetMaxId()+1];
    Maximum = new float [VariableNames->GetMaxId()+1];
    VectorLength = new int [VariableNames->GetMaxId()+1];
    this->MakeMeshFlag = 1;
    }

  output->DeepCopy(Mesh);  // If mesh has already been made copy it to output

  int first = 0;
  for (int j = 0; j <= VariableNames->GetMaxId(); j++)
    {
    if ( this->CellDataArraySelection->GetArraySetting(j) == 1 )
      {
      if (VariableComponents->GetValue(j) == 1)
        {
        GetVariableAtTimestep( j, this->TimeStep, CellDataArray[j]);
        }
      else
        {
        if ( !strcmp(CoordinateSystem,"CYLINDRICAL" ))
          {
          ConvertVectorFromCylindricalToCartesian( j-3, j-1);
          }
        FillVectorVariable( j-3, j-2, j-1, CellDataArray[j]);
        }
      if (first == 0)
        {
        output->GetCellData()->SetScalars(CellDataArray[j]);
        }
      else
        {
        output->GetCellData()->AddArray(CellDataArray[j]);
        }

      double tempRange[2];
      CellDataArray[j]->GetRange(tempRange, -1);
      Minimum[j] = tempRange[0];
      Maximum[j] = tempRange[1];
      VectorLength[j] = 1;
      first = 1;
      }
    }
}


//----------------------------------------------------------------------------
int vtkMFIXReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  if ( this->RequestInformationFlag == 0)
    {
    if ( !this->FileName )
      {
      this->NumberOfPoints = 0;
      this->NumberOfCells = 0;
      vtkErrorMacro("No filename specified");
      return 0;
      }

    SetProjectName(this->FileName);
    ReadRestartFile();
    CreateVariableNames();
    GetTimeSteps();
    CalculateMaxTimeStep();
    MakeTimeStepTable(VariableNames->GetMaxId()+1);
    GetNumberOfVariablesInSPXFiles();
    MakeSPXTimeStepIndexTable(VariableNames->GetMaxId()+1);

    for (int j = 0; j <= VariableNames->GetMaxId(); j++)
      {
      this->CellDataArraySelection->AddArray( VariableNames->GetValue(j));
      }

    this->NumberOfPoints = (IMaximum2+1)*(JMaximum2+1)*(KMaximum2+1);
    this->NumberOfCells = IJKMaximum2;
    this->NumberOfCellFields = VariableNames->GetMaxId()+1;
    this->NumberOfTimeSteps = this->MaximumTimestep;
    this->TimeStepRange[0] = 0;  
    this->TimeStepRange[1] = this->NumberOfTimeSteps-1;
    this->RequestInformationFlag = 1;
    this->GetAllTimes(outputVector);
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
void vtkMFIXReader::GetCellDataRange(int cellComp, int index, 
     float *min, float *max)
{
  if (index >= this->VectorLength[cellComp] || index < 0)
    {
    index = 0;  // if wrong index, set it to zero
    }
  *min = this->Minimum[cellComp];
  *max = this->Maximum[cellComp];
}

//----------------------------------------------------------------------------
void vtkMFIXReader::SetProjectName (char *infile) {
  int len = strlen(infile);
  strncpy(RunName, infile, len-4);
}

//----------------------------------------------------------------------------
void vtkMFIXReader::RestartVersionNumber(char* buffer)
{
  char s1[120];
  char s2[120];
  sscanf(buffer,"%s %s %f", s1, s2, &VersionNumber);
  strncpy(Version, buffer, 100);
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetInt(istream& in, int &val)
{
  in.read( (char*)&val,sizeof(int));
  SwapInt(val);
}

//----------------------------------------------------------------------------
void vtkMFIXReader::SwapInt(int &value)
{
  static char Swapped[4];
  int * Addr = &value;
  Swapped[0]=*((char*)Addr+3);
  Swapped[1]=*((char*)Addr+2);
  Swapped[2]=*((char*)Addr+1);
  Swapped[3]=*((char*)Addr  );
  value = *(reinterpret_cast<int*>(Swapped));
}

//----------------------------------------------------------------------------
void vtkMFIXReader::SwapDouble(double &value)
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

//----------------------------------------------------------------------------
void vtkMFIXReader::SwapFloat(float &value)
{
  static char Swapped[4];
  float * Addr = &value;

  Swapped[0]=*((char*)Addr+3);
  Swapped[1]=*((char*)Addr+2);
  Swapped[2]=*((char*)Addr+1);
  Swapped[3]=*((char*)Addr  );
  value = *(reinterpret_cast<float*>(Swapped));
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetDouble(istream& in, double& val)
{
  in.read( (char*)&val,sizeof(double));
  SwapDouble(val);
}

//----------------------------------------------------------------------------
void vtkMFIXReader::SkipBytes(istream& in, int n)
{
  in.read(DataBuffer,n); // maybe seekg instead
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetBlockOfDoubles(istream& in, vtkDoubleArray *v, int n)
{
  const int nr = 512/sizeof(double);
  double array[nr];
  int num_records;

  if ( n%nr == 0)
    {
    num_records = n/nr;
    }
  else
    {
    num_records = 1 + n/nr;
    }

  int c = 0;
  for (int i=0; i<num_records; ++i)
    {
    in.read( (char*)&array , 512 );
    for (int j=0; j<nr; ++j)
      {
      if (c < n) 
        {
        double temp = array[j];
        SwapDouble(temp);
        v->InsertValue( c, temp);
        ++c;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetBlockOfInts(istream& in, vtkIntArray *v, int n)
{
  const int nr = 512/sizeof(int);
  int array[nr];
  int num_records;

  if ( n%nr == 0)
    {
    num_records = n/nr;
    }
  else
    {
    num_records = 1 + n/nr;
    }

  int c = 0;
  for (int i = 0; i < num_records; ++i)
    {
    in.read( (char*)&array , 512 );
    for (int j=0; j<nr; ++j)
      {
      if (c < n)
        {
        int temp = array[j];
        SwapInt(temp);
        v->InsertValue( c, temp);
        ++c;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetBlockOfFloats(istream& in, vtkFloatArray *v, int n)
{
  const int nr = 512/sizeof(float);
  float array[nr];
  int num_records;

  if (!in)
    {
    cout << "Error opening file" << endl;
    }

  if ( n%nr == 0)
    {
    num_records = n/nr;
    }
  else
    {
    num_records = 1 + n/nr;
    }

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
        SwapFloat(temp);
        if ( Flag->GetValue(c) < 10) 
          {
          v->InsertValue(cnt, temp);
          cnt++;
          }
        ++c;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::ReadRestartFile()
{
  int i,l,n,lc;
  int DIMENSION_USR = 5;

  ifstream in(this->FileName,ios::binary);
  if (!in)
    {
    cout << "could not open file" << endl;
    return;
    }

  DataBuffer[512] = '\0';

  // version : record 1
  memset(DataBuffer,0,513);
  in.read(DataBuffer,512);
  RestartVersionNumber(DataBuffer);

  // skip 2 linesline : records 2 and 3
  in.read(DataBuffer,512);
  in.read(DataBuffer,512);

  // IMinimum1 etc : record 4
  memset(DataBuffer,0,513);

  if (Version == "RES = 01.00")
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetDouble(in,DeltaTime);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);

    // 15 ints ... 4 doubles = 92 bytes
    SkipBytes(in,420);
    }
  else if (Version == "RES = 01.01" || Version == "RES = 01.02")
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetInt(in,DimensionIc);    GetInt(in,DimensionBc);
    GetDouble(in,DeltaTime);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);

    // 17 ints ... 4 doubles = 100 bytes
    SkipBytes(in,412);
    }
  else if(Version == "RES = 01.03")
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetInt(in,DimensionIc); GetInt(in,DimensionBc);
    GetDouble(in,DeltaTime);
    GetDouble(in,XMinimum);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);

    // 17 ints ... 5 doubles = 108 bytes
    SkipBytes(in,404);
    }
  else if(Version == "RES = 01.04")
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetInt(in,DimensionIc); GetInt(in,DimensionBc);  GetInt(in,DimensionC);
    GetDouble(in,DeltaTime);
    GetDouble(in,XMinimum);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);

    // 18 ints ... 5 doubles = 112 bytes
    SkipBytes(in,400);
    }
  else if(Version == "RES = 01.05")
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetInt(in,DimensionIc); GetInt(in,DimensionBc);  GetInt(in,DimensionC);
    GetInt(in,DimensionIs);
    GetDouble(in,DeltaTime);
    GetDouble(in,XMinimum);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);

    // 19 ints ... 5 doubles = 116 bytes
    SkipBytes(in,396);
    }
  else
    {
    GetInt(in,IMinimum1);  GetInt(in,JMinimum1);   GetInt(in,KMinimum1);
    GetInt(in,IMaximum);   GetInt(in,JMaximum);    GetInt(in,KMaximum);
    GetInt(in,IMaximum1);  GetInt(in,JMaximum1);   GetInt(in,KMaximum1);
    GetInt(in,IMaximum2);  GetInt(in,JMaximum2);   GetInt(in,KMaximum2);
    GetInt(in,IJMaximum2); GetInt(in,IJKMaximum2); GetInt(in,MMAX);
    GetInt(in,DimensionIc); GetInt(in,DimensionBc);  GetInt(in,DimensionC);
    GetInt(in,DimensionIs);
    GetDouble(in,DeltaTime);
    GetDouble(in,XMinimum);
    GetDouble(in,XLength);  GetDouble(in,YLength);  GetDouble(in,ZLength);
    GetDouble(in,Ce); GetDouble(in,Cf); GetDouble(in,Phi); GetDouble(in,PhiW);

    // 19 ints ... 9 doubles = 148 bytes
    SkipBytes(in,364);
    }

  const int nr = 512/sizeof(float);

  if ( IJKMaximum2%nr == 0)
    {
    SPXRecordsPerTimestep = IJKMaximum2/nr;
    }
  else
    {
    SPXRecordsPerTimestep = 1 + IJKMaximum2/nr;
    }

  // C , C_name and nmax

  NMax->Resize(MMAX+1);
  for (int lc=0; lc<MMAX+1; ++lc)
    {
    NMax->InsertValue(lc, 1);
    }

  C->Resize(DimensionC);

  if (VersionNumber > 1.04)
    {
    GetBlockOfDoubles (in, C, DimensionC);

    for (lc=0; lc<DimensionC; ++lc) 
      {
      in.read(DataBuffer,512);  // c_name[]
      }

    if (VersionNumber < 1.12)
      {
      GetBlockOfInts(in, NMax,MMAX+1);
      }
    else
      {
      // what is the diff between this and above ??? 
      for (lc=0; lc<MMAX+1; ++lc) 
        {
        int temp;
        GetInt(in,temp);
        NMax->InsertValue(lc, temp);
        }

      SkipBytes(in,512-(MMAX+1)*sizeof(int));
      }
    }

  Dx->Resize(IMaximum2);
  Dy->Resize(JMaximum2);
  Dz->Resize(KMaximum2);

  GetBlockOfDoubles(in, Dx,IMaximum2);
  GetBlockOfDoubles(in, Dy,JMaximum2);
  GetBlockOfDoubles(in, Dz,KMaximum2);

  // RunName etc.

  memset(Units,0,17);
  memset(CoordinateSystem,0,17);

  in.read(DataBuffer,120);      // run_name , description
  in.read(Units,16);        // Units
  in.read(DataBuffer,16);       // run_type
  in.read(CoordinateSystem,16);  // CoordinateSystem 

  SkipBytes(in,512-168);

  char tmp[17];

  memset(tmp,0,17);

  int ic = 0;
  for (i=0; i<17; ++i)
    {
    if (Units[i] != ' ') 
      {
      tmp[ic++] = Units[i];
      }
    }

  memset(tmp,0,17);

  ic = 0;
  for (i=0; i<17; ++i)
    {
    if (CoordinateSystem[i] != ' ')
      {
      tmp[ic++] = CoordinateSystem[i];
      }
    }
  strcpy(CoordinateSystem,tmp);

  if (VersionNumber >= 1.04)
    {
    TempD->Resize(NMax->GetValue(0));
    GetBlockOfDoubles(in, TempD, NMax->GetValue(0));             // MW_g
    for (i=0; i<MMAX; ++i)
      {
      in.read(DataBuffer,512);  // MW_s
      }
    }

  in.read(DataBuffer,512);  // D_p etc.

  // read in the "DimensionIc" variables (and ignore ... not used by ani_mfix)
  TempI->Resize(DimensionIc);
  TempD->Resize(DimensionIc);

  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_x_w
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_x_e
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_y_s
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_y_n
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_z_b
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_z_t

  GetBlockOfInts(in, TempI,DimensionIc);  // ic_i_w
  GetBlockOfInts(in, TempI,DimensionIc);  // ic_i_e
  GetBlockOfInts(in, TempI,DimensionIc);  // ic_j_s
  GetBlockOfInts(in, TempI,DimensionIc);  // ic_j_n
  GetBlockOfInts(in, TempI,DimensionIc);  // ic_k_b
  GetBlockOfInts(in, TempI,DimensionIc);  // ic_k_t

  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_ep_g
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_p_g
  GetBlockOfDoubles(in, TempD,DimensionIc);  // ic_t_g

  if (VersionNumber < 1.15)
    {
    GetBlockOfDoubles(in,TempD,DimensionIc);  // ic_t_s(1,1)
    GetBlockOfDoubles(in,TempD,DimensionIc);  // ic_t_s(1,2) or ic_tmp 
    }

  if (VersionNumber >= 1.04)
    {
    for (int i=0; i<NMax->GetValue(0); ++i)
      {
      GetBlockOfDoubles(in,TempD,DimensionIc); // ic_x_g
      }
    }

  GetBlockOfDoubles(in,TempD,DimensionIc); // ic_u_g
  GetBlockOfDoubles(in,TempD,DimensionIc); // ic_v_g
  GetBlockOfDoubles(in,TempD,DimensionIc); // ic_w_g

  for (lc=0; lc<MMAX; ++lc)
    {
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_rop_s
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_u_s
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_v_s
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_w_s

    if (VersionNumber >= 1.15)
      {
      GetBlockOfDoubles(in,TempD,DimensionIc); // ic_t_s
      }

    if (VersionNumber >= 1.04)
      {
      for (n=0; n<NMax->GetValue(lc+1); ++n)
        {
        GetBlockOfDoubles(in,TempD,DimensionIc); // ic_x_s
        }
      }
    }

  // read in the "DimensionBc" variables (and ignore ... not used by ani_mfix)
  TempI->Resize(DimensionBc);
  TempD->Resize(DimensionBc);

  GetBlockOfDoubles(in,TempD,DimensionBc); // bc_x_w
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc_x_e
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc y s
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc y n
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc z b
  GetBlockOfDoubles(in,TempD,DimensionBc);  // bc z t
  GetBlockOfInts(in,TempI,DimensionBc);  // bc i w
  GetBlockOfInts(in,TempI,DimensionBc); // bc i e
  GetBlockOfInts(in,TempI,DimensionBc); // bc j s
  GetBlockOfInts(in,TempI,DimensionBc); // bc j n
  GetBlockOfInts(in,TempI,DimensionBc); // bc k b
  GetBlockOfInts(in,TempI,DimensionBc); // bc k t
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc ep g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc p g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc t g

  if (VersionNumber < 1.15)
    {
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_t_s(1,1)
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_t_s(1,1) or bc_tmp
    }

  if (VersionNumber >= 1.04)
    {
    for (int i=0; i<NMax->GetValue(0); ++i)
      {
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc_x_g
      }
    }

  GetBlockOfDoubles(in,TempD,DimensionBc); // bc u g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc v g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc w g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc ro g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc_rop_g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc volflow g
  GetBlockOfDoubles(in,TempD,DimensionBc); // bc massflow g

  for (lc=0; lc<MMAX; ++lc)
    {
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc rop s
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc u s
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc v s

    if (VersionNumber >= 1.04)
      {
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc w s

      if (VersionNumber >= 1.15)
        {
        GetBlockOfDoubles(in,TempD,DimensionBc); // bc t s
        }
      for (n=0; n<NMax->GetValue(lc+1); ++n)
        {
        GetBlockOfDoubles(in,TempD,DimensionBc); // bc x s
        }
      }
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc volflow s
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc massflow s
    }

  if (Version == "RES = 01.00")
    {
    l = 10;
    }
  else
    {
    l = DimensionBc;
    }

  for (lc=0; lc<l; ++lc)
    {
    in.read(DataBuffer,512); // BC TYPE
    }

  Flag->Resize(IJKMaximum2);
  GetBlockOfInts(in, Flag,IJKMaximum2);

  // DimensionIs varibles (not needed by ani_mfix)
  TempI->Resize(DimensionIs);
  TempD->Resize(DimensionIs);

  if (VersionNumber >= 1.04)
    {
    GetBlockOfDoubles(in,TempD,DimensionIs); // is x w
    GetBlockOfDoubles(in,TempD,DimensionIs); // is x e
    GetBlockOfDoubles(in,TempD,DimensionIs); // is y s
    GetBlockOfDoubles(in,TempD,DimensionIs); // is y n
    GetBlockOfDoubles(in,TempD,DimensionIs); // is z b
    GetBlockOfDoubles(in,TempD,DimensionIs); // is z t
    GetBlockOfInts(in,TempI,DimensionIs); // is i w
    GetBlockOfInts(in,TempI,DimensionIs); // is i e
    GetBlockOfInts(in,TempI,DimensionIs); // is j s
    GetBlockOfInts(in,TempI,DimensionIs); // is j n
    GetBlockOfInts(in,TempI,DimensionIs); // is k b
    GetBlockOfInts(in,TempI,DimensionIs); // is k t
    GetBlockOfDoubles(in,TempD,DimensionIs);  // is_pc(1,1)
    GetBlockOfDoubles(in,TempD,DimensionIs);  // is_pc(1,2)

    if (VersionNumber >= 1.07)
      {
      for (l=0; l<MMAX; ++l) GetBlockOfDoubles(in,TempD,DimensionIs);//is_vel_s
      }

    for (lc=0; lc<DimensionIs; ++lc)
      {
      in.read(DataBuffer,512); // is_type
      }
    }

  if (VersionNumber >= 1.08)
    {
    in.read(DataBuffer,512);
    }

  if (VersionNumber >= 1.09)
    {
    in.read(DataBuffer,512);

    if (VersionNumber >= 1.5)
      {
      GetInt(in,NumberOfSPXFilesUsed);
      SkipBytes(in,508);
      }

    for (lc=0; lc< NumberOfSPXFilesUsed; ++lc)
      {
      in.read(DataBuffer,512); // spx_dt
      }

    for (lc=0; lc<MMAX+1; ++lc)
      {
      in.read(DataBuffer,512);    // species_eq
      }

    TempD->Resize(DIMENSION_USR);

    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr_dt
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr x w
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr x e
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr y s
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr y n
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr z b
    GetBlockOfDoubles(in,TempD,DIMENSION_USR); // usr z t

    for (lc=0; lc<DIMENSION_USR; ++lc)
      {
      in.read(DataBuffer,512);    // usr_ext etc.
      }

    TempD->Resize(DimensionIc);
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_p_star
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_l_scale
    for (lc=0; lc<DimensionIc; ++lc)
      {
      in.read(DataBuffer,512);    // ic_type
      }

    TempD->Resize(DimensionBc);
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_dt_0
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_jet_g0
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_dt_h
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_jet_gh
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_dt_l
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_jet_gl
    }

  if (VersionNumber >= 1.1)
    {
    in.read(DataBuffer,512);  // mu_gmax
    }

  if (VersionNumber >= 1.11)
    {
    in.read(DataBuffer,512);  // x_ex , model_b
    }

  if (VersionNumber >= 1.12)
    {
    in.read(DataBuffer,512);   // p_ref , etc.
    in.read(DataBuffer,512);   // leq_it , leq_method

    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_hw_g
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_uw_g
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_vw_g
    GetBlockOfDoubles(in,TempD,DimensionBc); // bc_ww_g

    for (lc=0; lc<MMAX; ++lc)
      {
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc_hw_s
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc_uw_s
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc_vw_s
      GetBlockOfDoubles(in,TempD,DimensionBc); // bc_ww_s
      }
    }

  if (VersionNumber >= 1.13)
    {
    in.read(DataBuffer,512);    // momentum_x_eq , etc.
    }

  if (VersionNumber >= 1.14)
    {
    in.read(DataBuffer,512);    // detect_small
    }

  if (VersionNumber >= 1.15)
    {
    in.read(DataBuffer,512);    // k_g0 , etc.

    TempD->Resize(DimensionIc);

    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_gama_rg
    GetBlockOfDoubles(in,TempD,DimensionIc); // ic_t_rg

    for (lc=0; lc<MMAX; ++lc)
      {
      GetBlockOfDoubles(in,TempD,DimensionIc); // ic_gama_rs
      GetBlockOfDoubles(in,TempD,DimensionIc); // ic_t_rs
      }
    }

  if (VersionNumber >= 1.2)
    {
    in.read(DataBuffer,512); // norm_g , norm_s
    }

  if (VersionNumber >= 1.3)
    {
    GetInt(in,NumberOfScalars);
    SkipBytes(in,sizeof(double)); // tol_resid_scalar

    int DIM_tmp;
    GetInt(in,DIM_tmp);
    SkipBytes(in,512-sizeof(double)-2*sizeof(int));

    TempI->Resize(DIM_tmp);
    GetBlockOfInts(in,TempI,DIM_tmp);  // Phase4Scalar;
    }

  if (VersionNumber >= 1.5)
    {
    GetInt(in,NumberOfReactionRates);
    SkipBytes(in,508);
    }

  if (VersionNumber >= 1.5999)
    {
    int tmp;
    GetInt(in,tmp);
    SkipBytes(in,508);

    if (tmp != 0)
      {
      BkEpsilon = true;
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::CreateVariableNames()
{

  char fname[256];
  int cnt = 0;

  for (int i=0; i<NumberOfSPXFilesUsed; ++i)
    {
    for(int k = 0; k < (int)sizeof(fname); k++)
      {
      fname[k]=0;
      }
    strncpy(fname, this->FileName, strlen(this->FileName)-4);

    if (i==0)
      {
      strcat(fname, ".SP1");
      }
    else if (i==1)
      {
      strcat(fname, ".SP2");
      }
    else if (i==2)
      {
      strcat(fname, ".SP3");
      }
    else if (i==3)
      {
      strcat(fname, ".SP4");
      }
    else if (i==4)
      {
      strcat(fname, ".SP5");
      }
    else if (i==5)
      {
      strcat(fname, ".SP6");
      }
    else if (i==6)
      {
      strcat(fname, ".SP7");
      }
    else if (i==7)
      {
      strcat(fname, ".SP8");
      }
    else if (i==8)
      {
      strcat(fname, ".SP9");
      }
    else if (i==9)
      {
      strcat(fname, ".SPA");
      }
    else
      {
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
          VariableNames->InsertValue(cnt++,"EP_g");
          VariableIndexToSPX->InsertValue(cnt-1, 1);
          VariableComponents->InsertValue(cnt-1, 1);
          break;
          }

        case 2:
          {
          VariableNames->InsertValue(cnt++,"P_g");
          VariableIndexToSPX->InsertValue(cnt-1, 2);
          VariableComponents->InsertValue(cnt-1, 1);
          VariableNames->InsertValue(cnt++,"P_star");
          VariableIndexToSPX->InsertValue(cnt-1, 2);
          VariableComponents->InsertValue(cnt-1, 1);
          break;
          }

        case 3:
          {
          VariableNames->InsertValue(cnt++,"U_g");
          VariableIndexToSPX->InsertValue(cnt-1, 3);
          VariableComponents->InsertValue(cnt-1, 1);
          VariableNames->InsertValue(cnt++,"V_g");
          VariableIndexToSPX->InsertValue(cnt-1, 3);
          VariableComponents->InsertValue(cnt-1, 1);
          VariableNames->InsertValue(cnt++,"W_g");
          VariableIndexToSPX->InsertValue(cnt-1, 3);
          VariableComponents->InsertValue(cnt-1, 1);
          VariableNames->InsertValue(cnt++,"Gas Velocity");
          VariableIndexToSPX->InsertValue(cnt-1, 3);
          VariableComponents->InsertValue(cnt-1, 3);
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
            for(int k=0;k<(int)sizeof(us);k++)
              {
              us[k]=0;
              }
            for(int k=0;k<(int)sizeof(vs);k++)
              {
              vs[k]=0;
              }
            for(int k=0;k<(int)sizeof(ws);k++)
              {
              ws[k]=0;
              }
            for(int k=0;k<(int)sizeof(sv);k++)
              {
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
            VariableNames->InsertValue(cnt++, us);
            VariableIndexToSPX->InsertValue(cnt-1, 4);
            VariableComponents->InsertValue(cnt-1, 1);

            VariableNames->InsertValue(cnt++, vs);
            VariableIndexToSPX->InsertValue(cnt-1, 4);
            VariableComponents->InsertValue(cnt-1, 1);

            VariableNames->InsertValue(cnt++, ws);
            VariableIndexToSPX->InsertValue(cnt-1, 4);
            VariableComponents->InsertValue(cnt-1, 1);

            VariableNames->InsertValue(cnt++, sv);
            VariableIndexToSPX->InsertValue(cnt-1, 4);
            VariableComponents->InsertValue(cnt-1, 3);
            }
          break;
          }

        case 5:
          {
          char rops[120];
          char temp[120];

          for (int i=0; i<MMAX; ++i)
            {
            for(int k=0;k<(int)sizeof(rops);k++)
              {
              rops[k]=0;
              }
            strcpy(rops, "ROP_s_");
            sprintf(temp, "%d", i+1);
            strcat(rops, temp);
            VariableNames->InsertValue(cnt++, rops);
            VariableIndexToSPX->InsertValue(cnt-1, 5);
            VariableComponents->InsertValue(cnt-1, 1);
            }
          break;
          }

        case 6:
          {
          VariableNames->InsertValue(cnt++, "T_g");
          VariableIndexToSPX->InsertValue(cnt-1, 6);
          VariableComponents->InsertValue(cnt-1, 1);

          if (VersionNumber <= 1.15)
            {
            VariableNames->InsertValue(cnt++, "T_s_1");
            VariableIndexToSPX->InsertValue(cnt-1, 6);
            VariableComponents->InsertValue(cnt-1, 1);

            if (MMAX > 1)
              {
              VariableNames->InsertValue(cnt++, "T_s_2");
              VariableIndexToSPX->InsertValue(cnt-1, 6);
              VariableComponents->InsertValue(cnt-1, 1);
              }
            else
              {
              VariableNames->InsertValue(cnt++, "T_s_2_not_used");
              VariableIndexToSPX->InsertValue(cnt-1, 6);
              VariableComponents->InsertValue(cnt-1, 1);
              }
            }
          else
            {
            char ts[120];
            char temp[120];

            for (int i=0; i<MMAX; ++i)
              {
              for(int k=0;k<(int)sizeof(ts);k++)
                {
                ts[k]=0;
                }
              strcpy(ts, "T_s_");
              sprintf(temp, "%d", i+1);
              strcat(ts, temp);
              VariableNames->InsertValue(cnt++, ts);
              VariableIndexToSPX->InsertValue(cnt-1, 6);
              VariableComponents->InsertValue(cnt-1, 1);
              }
            }
          break;
          }

        case 7:
          {
          char var[120];
          char temp[120];

          for (int i=0; i<NMax->GetValue(0); ++i)
            {
            for (int k=0;k<(int)sizeof(var);k++)
              {
              var[k]=0;
              }
            strcpy(var, "X_g_");
            sprintf(temp, "%d", i+1);
            strcat(var, temp);
            VariableNames->InsertValue(cnt++, var);
            VariableIndexToSPX->InsertValue(cnt-1, 7);
            VariableComponents->InsertValue(cnt-1, 1);
            }

          for (int m=1; m<=MMAX; ++m)
            {
            for (int i=0; i<NMax->GetValue(m); ++i)
              {
              char temp1[120];
              char temp2[120];
              for (int k=0;k<(int)sizeof(var);k++)
                {
                var[k]=0;
                }
              strcpy(var, "X_s_");
              sprintf(temp1, "%d", m);
              sprintf(temp2, "%d", i+1);
              strcat(var, temp1);
              strcat(var, "_");
              strcat(var, temp2);
              VariableNames->InsertValue(cnt++, var);
              VariableIndexToSPX->InsertValue(cnt-1, 7);
              VariableComponents->InsertValue(cnt-1, 1);
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
            for (int k=0;k<(int)sizeof(var);k++)
              {
              var[k]=0;
              }
            strcpy(var, "Theta_m_");
            sprintf(temp, "%d", i+1);
            strcat(var, temp);
            VariableNames->InsertValue(cnt++, var);
            VariableIndexToSPX->InsertValue(cnt-1, 8);
            VariableComponents->InsertValue(cnt-1, 1);
            }
          break;
          }

        case 9:
          {
          char var[120];
          char temp[120];

          for (int i=0; i<NumberOfScalars; ++i)
            {
            for(int k=0;k<(int)sizeof(var);k++)
              {
              var[k]=0;
              }
            strcpy(var, "Scalar_");
            sprintf(temp, "%d", i+1);
            strcat(var, temp);
            VariableNames->InsertValue(cnt++, var);
            VariableIndexToSPX->InsertValue(cnt-1, 9);
            VariableComponents->InsertValue(cnt-1, 1);
            }
          break;
          }

        case 10:
          {
          char var[120];
          char temp[120];

          for (int i=0; i<NumberOfReactionRates; ++i)
            {
            for (int k=0;k<(int)sizeof(var);k++)
              {
              var[k]=0;
              }
            strcpy(var, "RRates_");
            sprintf(temp, "%d", i+1);
            strcat(var, temp);
            VariableNames->InsertValue(cnt++, var);
            VariableIndexToSPX->InsertValue(cnt-1, 10);
            VariableComponents->InsertValue(cnt-1, 1);
            }
          break;
          }

        case 11:
          {
          if (BkEpsilon)
            {
            VariableNames->InsertValue(cnt++, "k_turb_g");
            VariableIndexToSPX->InsertValue(cnt-1, 11);
            VariableComponents->InsertValue(cnt-1, 1);
            VariableNames->InsertValue(cnt++, "e_turb_g");
            VariableIndexToSPX->InsertValue(cnt-1, 11);
            VariableComponents->InsertValue(cnt-1, 1);
            }
          break;
          }

        default:
          {
          cout << "unknown SPx file : " << i << "\n";
          break;
          }
        }
      }
    else 
      {
      this->SpxFileExists->InsertValue(i, 0);
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetTimeSteps()
{
  int next_rec, num_rec;
  char fname[256];
  int cnt = 0;

  for (int i=0; i<NumberOfSPXFilesUsed; ++i)
    {
    for (int k=0;k<(int)sizeof(fname);k++)
      {
      fname[k]=0;
      }
    strncpy(fname, this->FileName, strlen(this->FileName)-4);
    if (i==0)
      {
      strcat(fname, ".SP1");
      }
    else if (i==1)
      {
      strcat(fname, ".SP2");
      }
    else if (i==2)
      {
      strcat(fname, ".SP3");
      }
    else if (i==3)
      {
      strcat(fname, ".SP4");
      }
    else if (i==4)
      {
      strcat(fname, ".SP5");
      }
    else if (i==5)
      {
      strcat(fname, ".SP6");
      }
    else if (i==6)
      {
      strcat(fname, ".SP7");
      }
    else if (i==7)
      {
      strcat(fname, ".SP8");
      }
    else if (i==8)
      {
      strcat(fname, ".SP9");
      }
    else if (i==9)
      {
      strcat(fname, ".SPA");
      }
    else
      {
      strcat(fname, ".SPB");
      }
    ifstream in(fname , ios::binary);

    int nvars=0;
    if (in) // file exists
      {
      in.clear();
      in.seekg( 1024, ios::beg ); 
      in.read( (char*)&next_rec,sizeof(int) );
      SwapInt(next_rec);
      in.read( (char*)&num_rec,sizeof(int) );
      SwapInt(num_rec);

      switch (i+1)
        {
        case 1: 
          {
          nvars = 1;
          break;
          }
        case 2:
          {
          nvars = 2;
          break;
          }
        case 3:
          {
          nvars = 4;
          break;
          }
        case 4:
          {
          nvars = 4*MMAX;
          break;
          }
        case 5:
          {
          nvars = MMAX;
          break;
          }
        case 6:
          {
          if (VersionNumber <= 1.15)
            {
            nvars = 3;
            }
          else
            {
            nvars = MMAX + 1;
            }
          break;
          }
        case 7:
          {
          nvars = NMax->GetValue(0);
          for (int m=0; m<MMAX; ++m)
            {
            nvars += NMax->GetValue(m);
            }
          break;
          }
        case 8:
          {
          nvars = MMAX;
          break;
          }
        case 9:
          {
          nvars = NumberOfScalars;
          break;
          }
        case 10:
          {
          nvars = NumberOfReactionRates;
          break;
          }
        case 11:
          {
          if (BkEpsilon)
            {
            nvars = 2;
            }
          break;
          }
        }

      for(int j=0; j<nvars; j++)
        {
        VariableTimesteps->InsertValue(cnt, (next_rec-4)/num_rec);
        cnt++;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::MakeTimeStepTable(int nvars)
{
  VariableTimestepTable->SetNumberOfComponents(nvars);

  for(int i=0; i<nvars; i++)
    {
    int ts_increment = MaximumTimestep/VariableTimesteps->GetValue(i);
    int ts = 1;
    for (int j=0; j<MaximumTimestep; j++)
      {
      VariableTimestepTable->InsertComponent(j, i, ts);
      ts_increment--;
      if (ts_increment <= 0)
        {
        ts_increment = MaximumTimestep/VariableTimesteps->GetValue(i);
        ts++;
        }
      if (ts > VariableTimesteps->GetValue(i)) 
        {
        ts = VariableTimesteps->GetValue(i);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetVariableAtTimestep(int vari , int tstep, 
  vtkFloatArray *v)
{
  // This routine opens and closes the file for each request.
  // Maybe keep all SPX files open, and just perform relative
  // moves to get to the correct location in the file
  // get filename that vaiable # vari is located in
  // assumptions : there are <10 solid phases,
  // <10 scalars and <10 ReactionRates (need to change this)

  char vname[256];
  strcpy(vname, VariableNames->GetValue(vari));
  int spx = VariableIndexToSPX->GetValue(vari);
  char fname[256];

  for(int k=0;k<(int)sizeof(fname);k++)
    {
    fname[k]=0;
    }

  strncpy(fname, this->FileName, strlen(this->FileName)-4);

  if (spx==1)
    {
    strcat(fname, ".SP1");
    }
  else if (spx==2)
    {
    strcat(fname, ".SP2");
    }
  else if (spx==3)
    {
    strcat(fname, ".SP3");
    }
  else if (spx==4)
    {
    strcat(fname, ".SP4");
    }
  else if (spx==5)
    {
    strcat(fname, ".SP5");
    }
  else if (spx==6)
    {
    strcat(fname, ".SP6");
    }
  else if (spx==7)
    {
    strcat(fname, ".SP7");
    }
  else if (spx==8)
    {
    strcat(fname, ".SP8");
    }
  else if (spx==9)
    {
    strcat(fname, ".SP9");
    }
  else if (spx==10)
    {
    strcat(fname, ".SPA");
    }
  else
    {
    strcat(fname, ".SPB");
    }

  int ind = (vari*MaximumTimestep) + tstep;
  int nBytesSkip = SPXTimestepIndexTable[ind];
  ifstream in(fname,ios::binary);
  in.seekg(nBytesSkip,ios::beg);
  GetBlockOfFloats (in, v, IJKMaximum2);
}

//----------------------------------------------------------------------------
void vtkMFIXReader::MakeSPXTimeStepIndexTable(int nvars)
{
  int SPXTimestepIndexTableSize = nvars * MaximumTimestep;
  SPXTimestepIndexTable = new int [SPXTimestepIndexTableSize];

  int timestep;
  int spx;
  int nvars_in_spx;

  for (int i=0; i<nvars; i++)
    {
    for (int j=0; j<MaximumTimestep; j++)
      {
      timestep = (int) VariableTimestepTable->GetComponent(j, i);
      spx = VariableIndexToSPX->GetValue(i);
      nvars_in_spx = SPXToNVarTable->GetValue(spx);
      int skip = VariableToSkipTable->GetValue(i);
      int index = (3*512) + (timestep-1) * 
        ((nvars_in_spx*SPXRecordsPerTimestep*512)+512) + 
        512 + (skip*SPXRecordsPerTimestep*512);
      int ind = (i*MaximumTimestep) + j;
      SPXTimestepIndexTable[ind] = index;
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::CalculateMaxTimeStep()
{
  this->MaximumTimestep = 0;
  for ( int i=0; i <= VariableNames->GetMaxId(); i++ )
    {
    if (VariableTimesteps->GetValue(i) > this->MaximumTimestep)
      {
      this->MaximumTimestep = VariableTimesteps->GetValue(i);
      }
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetNumberOfVariablesInSPXFiles()
{
  int spx_nvars = 0;
  int skip = 0;
  for (int j=1; j<NumberOfSPXFilesUsed; j++)
    {
    for(int i=0;i<VariableNames->GetMaxId()+1;i++)
      {
      if ((VariableIndexToSPX->GetValue(i) == j) 
        && (VariableComponents->GetValue(i) == 1))
        {
        spx_nvars++;
        VariableToSkipTable->InsertValue(i,skip);
        skip++;
        }
      }
    SPXToNVarTable->InsertValue(j, spx_nvars);
    spx_nvars = 0;
    skip = 0;
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::FillVectorVariable( int xindex, int yindex, 
  int zindex, vtkFloatArray *v)
{
  int range = CellDataArray[xindex]->GetMaxId();

  for(int i=0;i<=range;i++)
    {
    v->InsertComponent(i, 0, CellDataArray[xindex]->GetValue(i));
    v->InsertComponent(i, 1, CellDataArray[yindex]->GetValue(i));
    v->InsertComponent(i, 2, CellDataArray[zindex]->GetValue(i));
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::ConvertVectorFromCylindricalToCartesian( int xindex, 
  int zindex)
{
  int count = 0;
  float radius = 0.0;
  float y = 0.0;
  float theta = 0.0;
  int cnt=0;

  for (int k=0; k< KMaximum2; k++)
    {
    for (int j=0; j< JMaximum2; j++)
      {
      for (int i=0; i< IMaximum2; i++)
        {
        if ( Flag->GetValue(cnt) < 10 )
          {
          float ucart = (CellDataArray[xindex]->GetValue(count)*cos(theta)) -
            (CellDataArray[zindex]->GetValue(count)*sin(theta));
          float wcart = (CellDataArray[xindex]->GetValue(count)*sin(theta)) +
            (CellDataArray[zindex]->GetValue(count)*cos(theta));
          CellDataArray[xindex]->InsertValue(count, ucart);
          CellDataArray[zindex]->InsertValue(count, wcart);
          count++;
          }
        cnt++;
        radius = radius + Dx->GetValue(i);
        }
      radius = 0.0;
      y = y + Dy->GetValue(j);
      }
    y = 0.0;
    theta = theta + Dz->GetValue(k);
    }
}

//----------------------------------------------------------------------------
void vtkMFIXReader::GetAllTimes(vtkInformationVector *outputVector) 
{
  int max = 0;
  int maxvar = 0;

  for(int j=0; j<=VariableNames->GetMaxId(); j++)
    {
    int n = VariableTimesteps->GetValue(j);
    if (n > max)
      {
      max = n;
      maxvar = j;
      }
    }

  char fname[256];

  for(int k=0;k<(int)sizeof(fname);k++)
    {
    fname[k]=0;
    }
  strncpy(fname, this->FileName, strlen(this->FileName)-4);

  if (maxvar==0)
    {
    strcat(fname, ".SP1");
    }
  else if (maxvar==1)
    {
    strcat(fname, ".SP2");
    }
  else if (maxvar==2)
    {
    strcat(fname, ".SP3");
    }
  else if (maxvar==3)
    {
    strcat(fname, ".SP4");
    }
  else if (maxvar==4)
    {
    strcat(fname, ".SP5");
    }
  else if (maxvar==5)
    {
    strcat(fname, ".SP6");
    }
  else if (maxvar==6)
    {
    strcat(fname, ".SP7");
    }
  else if (maxvar==7)
    {
    strcat(fname, ".SP8");
    }
  else if (maxvar==8)
    {
    strcat(fname, ".SP9");
    }
  else if (maxvar==9)
    {
    strcat(fname, ".SPA");
    }
  else
    {
    strcat(fname, ".SPB");
    }

  ifstream tfile(fname , ios::binary);
  int numberOfVariablesInSPX = 
    SPXToNVarTable->GetValue(VariableIndexToSPX->GetValue(maxvar));
  int offset = 512-(int)sizeof(float) + 
    512*(numberOfVariablesInSPX*SPXRecordsPerTimestep);
  tfile.clear();
  tfile.seekg( 3*512, ios::beg ); // first time
  float time;
  double* steps = new double[this->NumberOfTimeSteps];

  for (int i = 0; i < this->NumberOfTimeSteps; i++)
    {
    tfile.read( (char*)&time,sizeof(float) );
    SwapFloat(time);
    steps[i] = (double)time;
    tfile.seekg(offset,ios::cur);
    }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
    steps, this->NumberOfTimeSteps);
}
