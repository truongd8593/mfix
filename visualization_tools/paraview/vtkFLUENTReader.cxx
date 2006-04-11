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
// Thanks to Brian W. Dotson (Department of Energy, National Energy 
//       Technology Laboratory) who developed this class.
//
// Please address all comments to Brian Dotson (Brian.Dotson@netl.doe.gov)

#include "vtkFLUENTReader.h"
#include "vtkDataArraySelection.h"
#include "vtkErrorCode.h"
#include "vtkUnstructuredGrid.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkByteSwap.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkByteSwap.h"
#include "vtkCellArray.h"
#include "vtkHexahedron.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkQuad.h"
#include "vtkTetra.h"
#include "vtkWedge.h"
#include "vtkPyramid.h"

vtkCxxRevisionMacro(vtkFLUENTReader, "$Revision$");
vtkStandardNewMacro(vtkFLUENTReader);

//----------------------------------------------------------------------------
vtkFLUENTReader::vtkFLUENTReader()
{
  this->FileName  = NULL;
  CreateVTKObjects();
}

//----------------------------------------------------------------------------
vtkFLUENTReader::~vtkFLUENTReader()
{
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->ReadFile(output);
  return 1;
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << endl;

  os << indent << "Number Of Cells: " << this->NumberOfCells << endl;
  os << indent << "Number Of Cell Fields: " 
     << this->NumberOfCellFields << endl;

  os << indent << "Number Of Cell Components: " 
     << this->NumberOfCellComponents << endl;

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::ReadFile(vtkUnstructuredGrid *output)
{
  output->Allocate();
  output->ShallowCopy(this->Mesh);
  this->Mesh->Delete();

  for ( int i=0; i < this->NumberOfVariables ; i++ ) {
    this->CellData[ i ]->Delete();
  }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  if(this->ObjectsFlag == 0){
    this->CreateVTKObjects();
  }

  if(!this->OpenCaseAndDataFiles()) {
    return 0;
  }

  this->ParseCaseFile();
  this->MakeFaceTreeParentTable();
  this->MakeCellTreeParentTable();
  this->LoadFaceKidFlags();
  this->LoadFaceParentFlags();
  this->LoadInterfaceFaceChildFlags();
  this->LoadNCGFaceChildFlags();
  this->LoadCellNumberOfFaces();
  this->LoadCellFaces();
  this->RemoveExtraFaces();
  this->LoadCellParentFlags();
  this->BuildCells();
  this->DataPass = 1;
  this->ParseDataFile();
  this->InitializeVariableNames();
  this->CellData = new vtkDoubleArray * [NumberOfVariables];

  for ( int i=0; i < this->NumberOfVariables ; i++ ) {
    int variableId = this->VariableIds->GetValue(i);
    int numberOfComponents = this->VariableSizes->GetValue(i);
    this->CellData[ i ] = vtkDoubleArray::New();
    this->CellData[ i ]->SetName(VariableNames[variableId]);
    this->CellData[ i ]->SetNumberOfComponents(numberOfComponents);
  }

  this->DataPass = 2;
  this->ParseDataFile();  // Getting Data

  int first = 0;
  for (int i=0; i<this->NumberOfVariables; i++ )
    {
    if((this->CellData[i]->GetNumberOfTuples() == this->NumberOfCells)
      && (this->CellData[i]->GetNumberOfComponents() < 6))
      {
      if(first == 0)
        {
        this->Mesh->GetCellData()->SetScalars(this->CellData[i]);
        } 
      else
        {
        this->Mesh->GetCellData()->AddArray(this->CellData[i]);
        }
      this->CellDataArraySelection->AddArray(this->CellData[ i ]->GetName());
      first = 1;
      this->NumberOfCellFields++;
      }
    }
  this->NumberOfCellArrays = this->NumberOfCellFields;
  this->Mesh->SetPoints(this->Points);
  this->DeleteVTKObjects();
  return 1;
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::OpenCaseAndDataFiles( void )
{
  int len = strlen(this->FileName);
  len = len -4;
  this->DataFileName = new char [256];
  strncpy(this->DataFileName, this->FileName, len);
  this->DataFileName[len] = '\0';
  strcat(this->DataFileName, ".dat");

  this->FileStream = new ifstream(this->FileName, ios::binary);
  this->DataFileStream = new ifstream(this->DataFileName, ios::binary);

  if (this->FileStream->fail()){
    cout << "Could Not Open Case File = " << this->FileName << endl;
    return(0);
  }

  if (this->DataFileStream->fail()){
    cout << "Could Not Open Data File = " << this->DataFileName << endl;
    return(0);
  }

  this->FileStream->seekg(0, ios::end); 
  this->CaseFileBufferLength = this->FileStream->tellg();
  this->FileStream->seekg(0, ios::beg);
  this->CaseFileBuffer = new char[this->CaseFileBufferLength];
  this->FileStream->read(this->CaseFileBuffer, this->CaseFileBufferLength);
  this->FileStream->close();

  this->DataFileStream->seekg(0, ios::end);
  this->DataFileBufferLength = this->DataFileStream->tellg();
  this->DataFileStream->seekg(0, ios::beg);
  this->DataFileBuffer = new char[this->DataFileBufferLength];
  this->DataFileStream->read(this->DataFileBuffer, this->DataFileBufferLength);
  this->DataFileStream->close();

  return(1);
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetCellDataRange( int cellComp, 
                                        int index, 
                                        float *min, 
                                        float *max)
{
  if (index >= this->VectorLength[cellComp] || index < 0)
    {
    index = 0;  // if wrong index, set it to zero
    }
  *min = this->Minimum[cellComp];
  *max = this->Maximum[cellComp];
}

//----------------------------------------------------------------------------
const char* vtkFLUENTReader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::SetCellArrayStatus(const char* name, int status)
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
void vtkFLUENTReader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::ParseCaseFile(void)
{
  int bufferIndex = 0;
  while(bufferIndex < this->CaseFileBufferLength)
    {
    if(this->CaseFileBuffer[bufferIndex] == '(')
      {
      int taskIndex = this->GetCaseIndex(bufferIndex);
      bufferIndex = this->ExecuteCaseTask(taskIndex, bufferIndex);
      }
      bufferIndex++;
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::MakeFaceTreeParentTable(void)
{
  for(int i = 0; i < this->NumberOfFaceTrees; i++)
    {
    if(this->FaceTreeParentFaceId1->GetValue(i) > this->LastFaceTreeParent)
      {
      this->LastFaceTreeParent = this->FaceTreeParentFaceId1->GetValue(i);
      }
    }

  for(int i=0; i <= this->LastFaceTreeParent; i++)
    {
    this->FaceTreeParentTable->InsertValue(i, 0);
    }

  int index = 0;
  for(int i = 0; i < this->NumberOfFaceTrees; i++)
    {
    for(int j = this->FaceTreeParentFaceId0->GetValue(i); 
      j <= this->FaceTreeParentFaceId1->GetValue(i); j++)
      {
      this->FaceTreeParentTable->InsertValue(j, index);
      index++;
      }
    }	
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::MakeCellTreeParentTable(void)
{
  for(int i = 0; i < this->NumberOfCellTrees; i++)
    {
    if(this->CellTreeParentCellId1->GetValue(i) > this->LastCellTreeParent)
      {
      this->LastCellTreeParent = this->CellTreeParentCellId1->GetValue(i);
      }
    }

  for(int i=0; i<=this->LastCellTreeParent; i++)
    {
    this->CellTreeParentTable->InsertValue(i, 0);
    }

  int index = 0;
  for(int i = 0; i < this->NumberOfCellTrees; i++)
    {
    for(int j = this->CellTreeParentCellId0->GetValue(i);
      j <= this->CellTreeParentCellId1->GetValue(i);j++)
      {
      this->CellTreeParentTable->InsertValue(j, index);
      index++;
      }
    }	
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadFaceKidFlags(void)
{
  // Initialize
  for(int i = 0; i <= this->NumberOfFaces; i++)
    {
    this->FaceKidFlags->InsertValue( i, 0);
    }

  for(int i=0;i<NumberOfFaceTrees;i++)
    {
    int startFace = this->FaceTreeParentFaceId0->GetValue(i);
    int endFace = this->FaceTreeParentFaceId1->GetValue(i);
    for(int j = startFace; j <= endFace; j++)
      {
      int startKid = this->FaceTreesKidsIndex->GetValue(
        this->FaceTreeParentTable->GetValue(j));
      int endKid = this->FaceTreesKidsIndex->GetValue(
        this->FaceTreeParentTable->GetValue(j))
        + this->FaceTreesNumberOfKids->GetValue(
        this->FaceTreeParentTable->GetValue(j));

        for(int k = startKid; k < endKid; k++)
        {
        int kid = this->FaceTreesKids->GetValue(k);
        this->FaceKidFlags->InsertValue( kid, 1);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadFaceParentFlags(void)
{
  // Initialize
  for(int i = 0; i <= this->NumberOfFaces; i++)
    {
    this->FaceParentFlags->InsertValue( i, 0);
    }

  for(int i = 0; i < this->NumberOfFaceTrees; i++)
    {
    int startFace = this->FaceTreeParentFaceId0->GetValue(i);
    int endFace = this->FaceTreeParentFaceId1->GetValue(i);
    for(int j = startFace; j <= endFace; j++)
      {
      this->FaceParentFlags->InsertValue( j, 1);
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadInterfaceFaceChildFlags(void)
{
  // Initialize Flag Array
  for(int i = 1; i <= this->NumberOfFaces; i++)
    {
    this->InterfaceFaceChildFlags->InsertValue(i,0);
    }

  for(int i = 0; i < this->NumberOfFaceParentChildren; i++)
    {
    int child = this->FaceParentsChildren->GetValue(i);
    this->InterfaceFaceChildFlags->InsertValue(child,1);
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadNCGFaceChildFlags(void)
{
  // Initialize Flag Array
  for(int i = 0; i <= this->NumberOfFaces; i++)
    {
    this->NCGFaceChildFlags->InsertValue(i,0);
    }

  for(int i = 0; i < this->NumberOfNCGFaces; i++)
    {
    int child = this->NCGFaceChild->GetValue(i);
    this->NCGFaceChildFlags->InsertValue(child,1);
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::BuildCells(void)
{
  int spinFace[6];
  for (int i = 0; i < 6; i++)
    {
    spinFace[i] = 0;
    }

  int node[8];
  for (int i = 0; i < 8; i++)
    {
    node[i] = 0;
    }

  int tempNode[30];
  for (int i = 0; i < 30; i++)
    {
    tempNode[i] = 0;
    }

  int face[6];

  for(int i=1;i<=NumberOfCells;i++)
    {
    for (int j = 0; j < 6; j++)
      {
      face[j] = (int) CellFacesClean->GetComponent(i, j);
      }

    for (int j = 0; j < 6; j++)
      {
      if ( (face[j]!=0) && ((int)FaceCells->GetComponent(face[j],0) == i))
        {
        spinFace[j] = 1;
        }
      else
        {
        spinFace[j] = -1;
        }
      }

    //*************************************
    //   Triangular Cell Type
    //*************************************

    if(CellTypes->GetValue(i) == 1)
      {
      int cnt = 0;
      for (int j = 0; j < 3; j++)
        {
        for (int k = 0; k < 2; k++)
          {
          tempNode[cnt] = (int)FaceNodes->GetComponent(face[j], k);
          }
        }

      if(spinFace[0] > 0)
        {
        node[0] = tempNode[0];
        node[1] = tempNode[1];
        }
      else
        {
        node[0] = tempNode[1];
        node[1] = tempNode[0];
        }

      if( (tempNode[2]!=node[0]) && (tempNode[2]!=node[1]) )
        {
        node[2] = tempNode[2];
        }
      else if ( (tempNode[3]!=node[0]) && (tempNode[3]!=node[1]) )
        {
        node[2] = tempNode[3];
        }
      else if ( (tempNode[4]!=node[0]) && (tempNode[4]!=node[1]) )
        {
        node[2] = tempNode[4];
        }
      else
        {
        node[2] = tempNode[5];
        }

      ATriangle->GetPointIds()->SetId( 0, node[0]);
      ATriangle->GetPointIds()->SetId( 1, node[1]);
      ATriangle->GetPointIds()->SetId( 2, node[2]);

      if(CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(ATriangle->GetCellType(), 
          ATriangle->GetPointIds());
        }
      }
    else if(CellTypes->GetValue(i) == 2)
      {
      //*************************************
      //   Tetrahedral Cell Type
      //*************************************
      tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
      tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
      tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);

      tempNode[3] = (int)FaceNodes->GetComponent(face[1], 0);
      tempNode[4] = (int)FaceNodes->GetComponent(face[1], 1);
      tempNode[5] = (int)FaceNodes->GetComponent(face[1], 2);

      tempNode[6] = (int)FaceNodes->GetComponent(face[2], 0);
      tempNode[7] = (int)FaceNodes->GetComponent(face[2], 1);
      tempNode[8] = (int)FaceNodes->GetComponent(face[2], 2);

      tempNode[9]  = (int)FaceNodes->GetComponent(face[3], 0);
      tempNode[10] = (int)FaceNodes->GetComponent(face[3], 1);
      tempNode[11] = (int)FaceNodes->GetComponent(face[3], 2);

      if (spinFace[0] > 0)
        {
        node[0] = tempNode[0];
        node[1] = tempNode[1];
        node[2] = tempNode[2];
        }
      else 
        {
        node[0] = tempNode[2];
        node[1] = tempNode[1];
        node[2] = tempNode[0];
        }

      if ( (tempNode[3]!=node[0]) && (tempNode[3]!=node[1]) && (tempNode[3]!=node[2]) ) 
        {
        node[3] = tempNode[3];
        }
      else if ( (tempNode[4]!=node[0]) && (tempNode[4]!=node[1]) && (tempNode[4]!=node[2]) ) 
        {
        node[3] = tempNode[4];
        }
      else if ( (tempNode[5]!=node[0]) && (tempNode[5]!=node[1]) && (tempNode[5]!=node[2]) ) 
        {
        node[3] = tempNode[5];
        } 
      else if ( (tempNode[6]!=node[0]) && (tempNode[6]!=node[1]) && (tempNode[6]!=node[2]) )
        {
        node[3] = tempNode[6];
	} 
      else if ( (tempNode[7]!=node[0]) && (tempNode[7]!=node[1]) && (tempNode[7]!=node[2]) ) 
        {
        node[3] = tempNode[7];
	} 
      else if ( (tempNode[8]!=node[0]) && (tempNode[8]!=node[1]) && (tempNode[8]!=node[2]) ) 
        {
        node[3] = tempNode[8];
        }
      else if ( (tempNode[9]!=node[0]) && (tempNode[9]!=node[1]) && (tempNode[9]!=node[2]) )
        {
        node[3] = tempNode[9];
        }
      else if ( (tempNode[10]!=node[0]) && (tempNode[10]!=node[1]) && (tempNode[10]!=node[2]) ) 
        {
        node[3] = tempNode[10];
        }
      else
        {
        node[3] = tempNode[11];
        }

      ATetra->GetPointIds()->SetId( 0, node[0]);
      ATetra->GetPointIds()->SetId( 1, node[1]);
      ATetra->GetPointIds()->SetId( 2, node[2]);
      ATetra->GetPointIds()->SetId( 3, node[3]);

      if (CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(ATetra->GetCellType(), ATetra->GetPointIds());
        }
      } 
    else if (CellTypes->GetValue(i) == 3)
      {
      //*************************************
      //   Quadrilateral Cell Type
      //*************************************

      tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
      tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
      tempNode[2] = (int)FaceNodes->GetComponent(face[1], 0);
      tempNode[3] = (int)FaceNodes->GetComponent(face[1], 1);
      tempNode[4] = (int)FaceNodes->GetComponent(face[2], 0);
      tempNode[5] = (int)FaceNodes->GetComponent(face[2], 1);
      tempNode[6] = (int)FaceNodes->GetComponent(face[3], 0);
      tempNode[7] = (int)FaceNodes->GetComponent(face[3], 1);

      if (spinFace[0] > 0)
        {
        node[0] = tempNode[0];
        node[1] = tempNode[1];
        }
      else
        {
        node[0] = tempNode[1];
        node[1] = tempNode[0];
        }

      if ( (tempNode[2]!=node[0]) && (tempNode[2]!=node[1]) && (tempNode[3]!=node[0]) && (tempNode[3]!=node[1]) )
        {
        if (spinFace[1] > 0)
          {
          node[2] = tempNode[2];
          node[3] = tempNode[3];
          }
        else 
          {
          node[2] = tempNode[3];
          node[3] = tempNode[2];
          }
        }

      if ( (tempNode[4]!=node[0]) && (tempNode[4]!=node[1]) && (tempNode[5]!=node[0]) && (tempNode[5]!=node[1]) )
        {
        if (spinFace[2] > 0)
          {
          node[2] = tempNode[4];
          node[3] = tempNode[5];
          }
        else 
          {
          node[2] = tempNode[5];
          node[3] = tempNode[4];
          }
        }

      if ( (tempNode[6]!=node[0]) && (tempNode[6]!=node[1]) && (tempNode[7]!=node[0]) && (tempNode[7]!=node[1]) )
        {
        if (spinFace[3] > 0)
          {
          node[2] = tempNode[6];
          node[3] = tempNode[7];
          }
        else
          {
          node[2] = tempNode[7];
          node[3] = tempNode[6];
          }
        }

      AQuad->GetPointIds()->SetId( 0, node[0]);
      AQuad->GetPointIds()->SetId( 1, node[1]);
      AQuad->GetPointIds()->SetId( 2, node[2]);
      AQuad->GetPointIds()->SetId( 3, node[3]);

      if (CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(AQuad->GetCellType(), AQuad->GetPointIds());
        }
      }
    else if (CellTypes->GetValue(i) == 4)
      {
      //*************************************
      //   Hexahedral Cell Type
      //*************************************
      int RightFace = 0;
      int LeftFace = 0;
      int FrontFace = 0;
      int BackFace = 0;
      int TopFace = 0;

      tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
      tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
      tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);
      tempNode[3] = (int)FaceNodes->GetComponent(face[0], 3);

      tempNode[4] = (int)FaceNodes->GetComponent(face[1], 0);
      tempNode[5] = (int)FaceNodes->GetComponent(face[1], 1);
      tempNode[6] = (int)FaceNodes->GetComponent(face[1], 2);
      tempNode[7] = (int)FaceNodes->GetComponent(face[1], 3);

      tempNode[8] =  (int)FaceNodes->GetComponent(face[2], 0);
      tempNode[9] =  (int)FaceNodes->GetComponent(face[2], 1);
      tempNode[10] = (int)FaceNodes->GetComponent(face[2], 2);
      tempNode[11] = (int)FaceNodes->GetComponent(face[2], 3);

      tempNode[12] = (int)FaceNodes->GetComponent(face[3], 0);
      tempNode[13] = (int)FaceNodes->GetComponent(face[3], 1);
      tempNode[14] = (int)FaceNodes->GetComponent(face[3], 2);
      tempNode[15] = (int)FaceNodes->GetComponent(face[3], 3);

      tempNode[16] = (int)FaceNodes->GetComponent(face[4], 0);
      tempNode[17] = (int)FaceNodes->GetComponent(face[4], 1);
      tempNode[18] = (int)FaceNodes->GetComponent(face[4], 2);
      tempNode[19] = (int)FaceNodes->GetComponent(face[4], 3);

      tempNode[20] = (int)FaceNodes->GetComponent(face[5], 0);
      tempNode[21] = (int)FaceNodes->GetComponent(face[5], 1);
      tempNode[22] = (int)FaceNodes->GetComponent(face[5], 2);
      tempNode[23] = (int)FaceNodes->GetComponent(face[5], 3);

      if (spinFace[0] > 0)
        {
        node[0] = tempNode[0];
        node[1] = tempNode[1];
        node[2] = tempNode[2];
        node[3] = tempNode[3];
        }
      else
        {
        node[0] = tempNode[3];
        node[1] = tempNode[2];
        node[2] = tempNode[1];
        node[3] = tempNode[0];
        }

      int FF = 0;
      int TF = 0;
      int RF = 0;
      int LF = 0;
      int BF = 0;

      int FFN0 = 0;
      int FFN1 = 0;
      int FFN2 = 0;
      int FFN3 = 0;

      int TFN0 = 0;
      int TFN1 = 0;
      int TFN2 = 0;
      int TFN3 = 0;

      int RFN0 = 0;
      int RFN1 = 0;
      int RFN2 = 0;
      int RFN3 = 0;

      int LFN0 = 0;
      int LFN1 = 0;
      int LFN2 = 0;
      int LFN3 = 0;

      int BFN0 = 0;
      int BFN1 = 0;
      int BFN2 = 0;
      int BFN3 = 0;

      if (((tempNode[4]==node[0])||(tempNode[5]==node[0])||(tempNode[6]==node[0])||(tempNode[7]==node[0])) 
        && ((tempNode[4]==node[1])||(tempNode[5]==node[1])||(tempNode[6]==node[1])||(tempNode[7]==node[1])) )
        {
        RightFace = 1;
        RF = face[1];
        RFN0 = tempNode[4];
        RFN1 = tempNode[5];
        RFN2 = tempNode[6];
        RFN3 = tempNode[7];
        }
      else if (((tempNode[4]==node[0])||(tempNode[5]==node[0])||(tempNode[6]==node[0])||(tempNode[7]==node[0]))
        && ((tempNode[4]==node[3])||(tempNode[5]==node[3])||(tempNode[6]==node[3])||(tempNode[7]==node[3])) )
        {
        FrontFace = 1;
        FF = face[1];
        FFN0 = tempNode[4];
        FFN1 = tempNode[5];
        FFN2 = tempNode[6];
        FFN3 = tempNode[7];
        }
      else if (((tempNode[4]==node[2])||(tempNode[5]==node[2])||(tempNode[6]==node[2])||(tempNode[7]==node[2])) 
        && ((tempNode[4]==node[3])||(tempNode[5]==node[3])||(tempNode[6]==node[3])||(tempNode[7]==node[3])) )
        {
        LeftFace = 1;
        LF = face[1];
        LFN0 = tempNode[4];
        LFN1 = tempNode[5];
        LFN2 = tempNode[6];
        LFN3 = tempNode[7];
        }
      else if (((tempNode[4]==node[1])||(tempNode[5]==node[1])||(tempNode[6]==node[1])||(tempNode[7]==node[1])) 
        && ((tempNode[4]==node[2])||(tempNode[5]==node[2])||(tempNode[6]==node[2])||(tempNode[7]==node[2])) ) 
        {
        BackFace = 1;
        BF = face[1];
        BFN0 = tempNode[4];
        BFN1 = tempNode[5];
        BFN2 = tempNode[6];
        BFN3 = tempNode[7];
        }
      else
        {
        TopFace = 1;
        TF = face[1];
        TFN0 = tempNode[4];
        TFN1 = tempNode[5];
        TFN2 = tempNode[6];
        TFN3 = tempNode[7];
        }

      if (((tempNode[8]==node[0])||(tempNode[9]==node[0])||(tempNode[10]==node[0])||(tempNode[11]==node[0])) 
        && ((tempNode[8]==node[1])||(tempNode[9]==node[1])||(tempNode[10]==node[1])||(tempNode[11]==node[1])))
        {
        RightFace = 2;
        RF = face[2];
        RFN0 = tempNode[8];
        RFN1 = tempNode[9];
        RFN2 = tempNode[10];
        RFN3 = tempNode[11];
        }
      else if (((tempNode[8]==node[0])||(tempNode[9]==node[0])||(tempNode[10]==node[0])||(tempNode[11]==node[0])) 
        && ((tempNode[8]==node[3])||(tempNode[9]==node[3])||(tempNode[10]==node[3])||(tempNode[11]==node[3])))
        {
        FrontFace = 2;
        FF = face[2];
        FFN0 = tempNode[8];
        FFN1 = tempNode[9];
        FFN2 = tempNode[10];
        FFN3 = tempNode[11];
        }
      else if (((tempNode[8]==node[2])||(tempNode[9]==node[2])||(tempNode[10]==node[2])||(tempNode[11]==node[2]))
        && ((tempNode[8]==node[3])||(tempNode[9]==node[3])||(tempNode[10]==node[3])||(tempNode[11]==node[3])))
        {
        LeftFace = 2;
        LF = face[2];
        LFN0 = tempNode[8];
        LFN1 = tempNode[9];
        LFN2 = tempNode[10];
        LFN3 = tempNode[11];
        }
      else if (((tempNode[8]==node[1])||(tempNode[9]==node[1])||(tempNode[10]==node[1])||(tempNode[11]==node[1]))
        && ((tempNode[8]==node[2])||(tempNode[9]==node[2])||(tempNode[10]==node[2])||(tempNode[11]==node[2])))
        {
        BackFace = 2;
        BF = face[2];
        BFN0 = tempNode[8];
        BFN1 = tempNode[9];
        BFN2 = tempNode[10];
        BFN3 = tempNode[11];
        }
      else
        {
        TopFace = 2;
        TF = face[2];
        TFN0 = tempNode[8];
        TFN1 = tempNode[9];
        TFN2 = tempNode[10];
        TFN3 = tempNode[11];
        }

      if (((tempNode[12]==node[0])||(tempNode[13]==node[0])||(tempNode[14]==node[0])||(tempNode[15]==node[0]))
        && ((tempNode[12]==node[1])||(tempNode[13]==node[1])||(tempNode[14]==node[1])||(tempNode[15]==node[1])))
        {
        RightFace = 3;
        RF = face[3];
        RFN0 = tempNode[12];
        RFN1 = tempNode[13];
        RFN2 = tempNode[14];
        RFN3 = tempNode[15];
        }
      else if (((tempNode[12]==node[0])||(tempNode[13]==node[0])||(tempNode[14]==node[0])||(tempNode[15]==node[0]))
        && ((tempNode[12]==node[3])||(tempNode[13]==node[3])||(tempNode[14]==node[3])||(tempNode[15]==node[3])))
        {
        FrontFace = 3;
        FF = face[3];
        FFN0 = tempNode[12];
        FFN1 = tempNode[13];
        FFN2 = tempNode[14];
        FFN3 = tempNode[15];
        }
      else if (((tempNode[12]==node[2])||(tempNode[13]==node[2])||(tempNode[14]==node[2])||(tempNode[15]==node[2]))
        && ((tempNode[12]==node[3])||(tempNode[13]==node[3])||(tempNode[14]==node[3])||(tempNode[15]==node[3])))
        {
        LeftFace = 3;
        LF = face[3];
        LFN0 = tempNode[12];
        LFN1 = tempNode[13];
        LFN2 = tempNode[14];
        LFN3 = tempNode[15];
        }
      else if (((tempNode[12]==node[1])||(tempNode[13]==node[1])||(tempNode[14]==node[1])||(tempNode[15]==node[1]))
        && ((tempNode[12]==node[2])||(tempNode[13]==node[2])||(tempNode[14]==node[2])||(tempNode[15]==node[2])))
        {
        BackFace = 3;
        BF = face[3];
        BFN0 = tempNode[12];
        BFN1 = tempNode[13];
        BFN2 = tempNode[14];
        BFN3 = tempNode[15];
        }
      else
        {
        TopFace = 3;
        TF = face[3];
        TFN0 = tempNode[12];
        TFN1 = tempNode[13];
        TFN2 = tempNode[14];
        TFN3 = tempNode[15];
        }

      if (((tempNode[16]==node[0])||(tempNode[17]==node[0])||(tempNode[18]==node[0])||(tempNode[19]==node[0]))
        && ((tempNode[16]==node[1])||(tempNode[17]==node[1])||(tempNode[18]==node[1])||(tempNode[19]==node[1])))
        {
        RightFace = 4;
        RF = face[4];
        RFN0 = tempNode[16];
        RFN1 = tempNode[17];
        RFN2 = tempNode[18];
        RFN3 = tempNode[19];
        }
      else if (((tempNode[16]==node[0])||(tempNode[17]==node[0])||(tempNode[18]==node[0])||(tempNode[19]==node[0]))
        && ((tempNode[16]==node[3])||(tempNode[17]==node[3])||(tempNode[18]==node[3])||(tempNode[19]==node[3])))
        {
        FrontFace = 4;
        FF = face[4];
        FFN0 = tempNode[16];
        FFN1 = tempNode[17];
        FFN2 = tempNode[18];
        FFN3 = tempNode[19];
        }
      else if (((tempNode[16]==node[2])||(tempNode[17]==node[2])||(tempNode[18]==node[2])||(tempNode[19]==node[2]))
        && ((tempNode[16]==node[3])||(tempNode[17]==node[3])||(tempNode[18]==node[3])||(tempNode[19]==node[3])))
        {
        LeftFace = 4;
        LF = face[4];
        LFN0 = tempNode[16];
        LFN1 = tempNode[17];
        LFN2 = tempNode[18];
        LFN3 = tempNode[19];
        }
      else if (((tempNode[16]==node[1])||(tempNode[17]==node[1])||(tempNode[18]==node[1])||(tempNode[19]==node[1]))
        && ((tempNode[16]==node[2])||(tempNode[17]==node[2])||(tempNode[18]==node[2])||(tempNode[19]==node[2])))
        {
        BackFace = 4;
        BF = face[4];
        BFN0 = tempNode[16];
        BFN1 = tempNode[17];
        BFN2 = tempNode[18];
        BFN3 = tempNode[19];
        }
      else
        {
        TopFace = 4;
        TF = face[4];
        TFN0 = tempNode[16];
        TFN1 = tempNode[17];
        TFN2 = tempNode[18];
        TFN3 = tempNode[19];
        }

      if (((tempNode[20]==node[0])||(tempNode[21]==node[0])||(tempNode[22]==node[0])||(tempNode[23]==node[0]))
        && ((tempNode[20]==node[1])||(tempNode[21]==node[1])||(tempNode[22]==node[1])||(tempNode[23]==node[1])))
        {
        RightFace = 5;
        RF = face[5];
        RFN0 = tempNode[20];
        RFN1 = tempNode[21];
        RFN2 = tempNode[22];
        RFN3 = tempNode[23];
        }
      else if (((tempNode[20]==node[0])||(tempNode[21]==node[0])||(tempNode[22]==node[0])||(tempNode[23]==node[0]))
        && ((tempNode[20]==node[3])||(tempNode[21]==node[3])||(tempNode[22]==node[3])||(tempNode[23]==node[3])))
        {
        FrontFace = 5;
        FF = face[5];
        FFN0 = tempNode[20];
        FFN1 = tempNode[21];
        FFN2 = tempNode[22];
        FFN3 = tempNode[23];
        }
      else if (((tempNode[20]==node[2])||(tempNode[21]==node[2])||(tempNode[22]==node[2])||(tempNode[23]==node[2])) 
        && ((tempNode[20]==node[3])||(tempNode[21]==node[3])||(tempNode[22]==node[3])||(tempNode[23]==node[3])))
        {
        LeftFace = 5;
        LF = face[5];
        LFN0 = tempNode[20];
        LFN1 = tempNode[21];
        LFN2 = tempNode[22];
        LFN3 = tempNode[23];
        }
      else if (((tempNode[20]==node[1])||(tempNode[21]==node[1])||(tempNode[22]==node[1])||(tempNode[23]==node[1])) 
        && ((tempNode[20]==node[2])||(tempNode[21]==node[2])||(tempNode[22]==node[2])||(tempNode[23]==node[2])))
        {
        BackFace = 5;
        BF = face[5];
        BFN0 = tempNode[20];
        BFN1 = tempNode[21];
        BFN2 = tempNode[22];
        BFN3 = tempNode[23];
        }
      else
        {
        TopFace = 5;
        TF = face[5];
        TFN0 = tempNode[20];
        TFN1 = tempNode[21];
        TFN2 = tempNode[22];
        TFN3 = tempNode[23];
        }

      if (((TFN0==RFN0)||(TFN0==RFN1)||(TFN0==RFN2)||(TFN0==RFN3))
        &&((TFN0==FFN0)||(TFN0==FFN1)||(TFN0==FFN2)||(TFN0==FFN3)))
        {
        node[4] = TFN0;
        }
      else if (((TFN0==RFN0)||(TFN0==RFN1)||(TFN0==RFN2)||(TFN0==RFN3))
        &&((TFN0==BFN0)||(TFN0==BFN1)||(TFN0==BFN2)||(TFN0==BFN3)))
        {
        node[5] = TFN0;
        }
      else if (((TFN0==LFN0)||(TFN0==LFN1)||(TFN0==LFN2)||(TFN0==LFN3))
        &&((TFN0==BFN0)||(TFN0==BFN1)||(TFN0==BFN2)||(TFN0==BFN3)))
        {
        node[6] = TFN0;
        }
      else 
        {
        node[7] = TFN0;
        }

      if (((TFN1==RFN0)||(TFN1==RFN1)||(TFN1==RFN2)||(TFN1==RFN3))
        &&((TFN1==FFN0)||(TFN1==FFN1)||(TFN1==FFN2)||(TFN1==FFN3)))
        {
        node[4] = TFN1;
        }
      else if (((TFN1==RFN0)||(TFN1==RFN1)||(TFN1==RFN2)||(TFN1==RFN3))
        &&((TFN1==BFN0)||(TFN1==BFN1)||(TFN1==BFN2)||(TFN1==BFN3)))
        {
        node[5] = TFN1;
        }
      else if (((TFN1==LFN0)||(TFN1==LFN1)||(TFN1==LFN2)||(TFN1==LFN3))
        &&((TFN1==BFN0)||(TFN1==BFN1)||(TFN1==BFN2)||(TFN1==BFN3)))
        {
        node[6] = TFN1;
        }
      else
        {
        node[7] = TFN1;
        }

      if (((TFN2==RFN0)||(TFN2==RFN1)||(TFN2==RFN2)||(TFN2==RFN3))
        &&((TFN2==FFN0)||(TFN2==FFN1)||(TFN2==FFN2)||(TFN2==FFN3)))
        {
        node[4] = TFN2;
        }
      else if (((TFN2==RFN0)||(TFN2==RFN1)||(TFN2==RFN2)||(TFN2==RFN3))
        &&((TFN2==BFN0)||(TFN2==BFN1)||(TFN2==BFN2)||(TFN2==BFN3)))
        {
        node[5] = TFN2;
        }
      else if (((TFN2==LFN0)||(TFN2==LFN1)||(TFN2==LFN2)||(TFN2==LFN3))
        &&((TFN2==BFN0)||(TFN2==BFN1)||(TFN2==BFN2)||(TFN2==BFN3)))
        {
        node[6] = TFN2;
        }
      else
        {
        node[7] = TFN2;
        }

      if (((TFN3==RFN0)||(TFN3==RFN1)||(TFN3==RFN2)||(TFN3==RFN3))
        &&((TFN3==FFN0)||(TFN3==FFN1)||(TFN3==FFN2)||(TFN3==FFN3)))
        {
        node[4] = TFN3;
        }
      else if (((TFN3==RFN0)||(TFN3==RFN1)||(TFN3==RFN2)||(TFN3==RFN3)) 
        &&((TFN3==BFN0)||(TFN3==BFN1)||(TFN3==BFN2)||(TFN3==BFN3)))
        {
        node[5] = TFN3;
        }
      else if (((TFN3==LFN0)||(TFN3==LFN1)||(TFN3==LFN2)||(TFN3==LFN3))
        &&((TFN3==BFN0)||(TFN3==BFN1)||(TFN3==BFN2)||(TFN3==BFN3)))
        {
        node[6] = TFN3;
        }
      else
        {
        node[7] = TFN3;
        }

      AHexahedron->GetPointIds()->SetId( 0, node[0]);
      AHexahedron->GetPointIds()->SetId( 1, node[1]);
      AHexahedron->GetPointIds()->SetId( 2, node[2]);
      AHexahedron->GetPointIds()->SetId( 3, node[3]);
      AHexahedron->GetPointIds()->SetId( 4, node[4]);
      AHexahedron->GetPointIds()->SetId( 5, node[5]);
      AHexahedron->GetPointIds()->SetId( 6, node[6]);
      AHexahedron->GetPointIds()->SetId( 7, node[7]);

      if (CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(AHexahedron->GetCellType(), 
          AHexahedron->GetPointIds());
        }
      }
    else if (CellTypes->GetValue(i) == 5)
      {
      //*************************************
      //   Pyramid Cell Type
      //*************************************
      int BF;
      if (FaceTypes->GetValue(face[0]) == 4)
        {
        BF = face[0];
        if (spinFace[0] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[0], 0);
          node[1] = (int)FaceNodes->GetComponent(face[0], 1);
          node[2] = (int)FaceNodes->GetComponent(face[0], 2);
          node[3] = (int)FaceNodes->GetComponent(face[0], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[0], 0);
          node[2] = (int)FaceNodes->GetComponent(face[0], 1);
          node[1] = (int)FaceNodes->GetComponent(face[0], 2);
          node[0] = (int)FaceNodes->GetComponent(face[0], 3);
          }
        tempNode[0] = (int)FaceNodes->GetComponent(face[1], 0);
        tempNode[1] = (int)FaceNodes->GetComponent(face[1], 1);
        tempNode[2] = (int)FaceNodes->GetComponent(face[1], 2);
        }
      else if (FaceTypes->GetValue(face[1]) == 4)
        {
        BF = face[1];
        if (spinFace[1] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[1], 0);
          node[1] = (int)FaceNodes->GetComponent(face[1], 1);
          node[2] = (int)FaceNodes->GetComponent(face[1], 2);
          node[3] = (int)FaceNodes->GetComponent(face[1], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[1], 0);
          node[2] = (int)FaceNodes->GetComponent(face[1], 1);
          node[1] = (int)FaceNodes->GetComponent(face[1], 2);
          node[0] = (int)FaceNodes->GetComponent(face[1], 3);
          }
          tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
          tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
          tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);
          }
      else if (FaceTypes->GetValue(face[2]) == 4)
        {
        BF = face[2];
        if (spinFace[2] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[2], 0);
          node[1] = (int)FaceNodes->GetComponent(face[2], 1);
          node[2] = (int)FaceNodes->GetComponent(face[2], 2);
          node[3] = (int)FaceNodes->GetComponent(face[2], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[2], 0);
          node[2] = (int)FaceNodes->GetComponent(face[2], 1);
          node[1] = (int)FaceNodes->GetComponent(face[2], 2);
          node[0] = (int)FaceNodes->GetComponent(face[2], 3);
          }
        tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
        tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
        tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);
        }
      else if (FaceTypes->GetValue(face[3]) == 4)
        {
        BF = face[3];
        if (spinFace[3] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[3], 0);
          node[1] = (int)FaceNodes->GetComponent(face[3], 1);
          node[2] = (int)FaceNodes->GetComponent(face[3], 2);
          node[3] = (int)FaceNodes->GetComponent(face[3], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[3], 0);
          node[2] = (int)FaceNodes->GetComponent(face[3], 1);
          node[1] = (int)FaceNodes->GetComponent(face[3], 2);
          node[0] = (int)FaceNodes->GetComponent(face[3], 3);
          }
        tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
        tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
        tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);
        }
      else
        {
        BF = face[4];
        if (spinFace[4] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[4], 0);
          node[1] = (int)FaceNodes->GetComponent(face[4], 1);
          node[2] = (int)FaceNodes->GetComponent(face[4], 2);
          node[3] = (int)FaceNodes->GetComponent(face[4], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[4], 0);
          node[2] = (int)FaceNodes->GetComponent(face[4], 1);
          node[1] = (int)FaceNodes->GetComponent(face[4], 2);
          node[0] = (int)FaceNodes->GetComponent(face[4], 3);
          }
        tempNode[0] = (int)FaceNodes->GetComponent(face[0], 0);
        tempNode[1] = (int)FaceNodes->GetComponent(face[0], 1);
        tempNode[2] = (int)FaceNodes->GetComponent(face[0], 2);
        }

      if ((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])&&(tempNode[0]!=node[2])&&(tempNode[0]!=node[3]))
        {
        node[4] = tempNode[0];
        }
      else if ((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])&&(tempNode[1]!=node[2])&&(tempNode[1]!=node[3]))
        {
        node[4] = tempNode[1];
        }
      else
        {
        node[4] = tempNode[2];
        }

      APyramid->GetPointIds()->SetId( 0, node[0]);
      APyramid->GetPointIds()->SetId( 1, node[1]);
      APyramid->GetPointIds()->SetId( 2, node[2]);
      APyramid->GetPointIds()->SetId( 3, node[3]);
      APyramid->GetPointIds()->SetId( 4, node[4]);

      if (CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(APyramid->GetCellType(),
          APyramid->GetPointIds());
        }

      }
    else if (CellTypes->GetValue(i) == 6)
      {
      //*************************************
      //   Wedge Cell Type
      //*************************************
      int BF;

      if (FaceTypes->GetValue(face[0]) == 4)
        {
        BF = face[0];
        if (spinFace[0] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[0], 0);
          node[1] = (int)FaceNodes->GetComponent(face[0], 1);
          node[4] = (int)FaceNodes->GetComponent(face[0], 2);
          node[3] = (int)FaceNodes->GetComponent(face[0], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[0], 0);
          node[4] = (int)FaceNodes->GetComponent(face[0], 1);
          node[1] = (int)FaceNodes->GetComponent(face[0], 2);
          node[0] = (int)FaceNodes->GetComponent(face[0], 3);
          }
        }
      else if (FaceTypes->GetValue(face[1]) == 4)
        {
        BF = face[1];
        if (spinFace[1] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[1], 0);
          node[1] = (int)FaceNodes->GetComponent(face[1], 1);
          node[4] = (int)FaceNodes->GetComponent(face[1], 2);
          node[3] = (int)FaceNodes->GetComponent(face[1], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[1], 0);
          node[4] = (int)FaceNodes->GetComponent(face[1], 1);
          node[1] = (int)FaceNodes->GetComponent(face[1], 2);
          node[0] = (int)FaceNodes->GetComponent(face[1], 3);
          }
        }
      else if  (FaceTypes->GetValue(face[2]) == 4)
        {
        BF = face[2];
        if (spinFace[2] > 0) 
          {
          node[0] = (int)FaceNodes->GetComponent(face[2], 0);
          node[1] = (int)FaceNodes->GetComponent(face[2], 1);
          node[4] = (int)FaceNodes->GetComponent(face[2], 2);
          node[3] = (int)FaceNodes->GetComponent(face[2], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[2], 0);
          node[4] = (int)FaceNodes->GetComponent(face[2], 1);
          node[1] = (int)FaceNodes->GetComponent(face[2], 2);
          node[0] = (int)FaceNodes->GetComponent(face[2], 3);
          }
        }
      else if (FaceTypes->GetValue(face[3]) == 4)
        {
        BF = face[3];
        if (spinFace[3] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[3], 0);
          node[1] = (int)FaceNodes->GetComponent(face[3], 1);
          node[4] = (int)FaceNodes->GetComponent(face[3], 2);
          node[3] = (int)FaceNodes->GetComponent(face[3], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[3], 0);
          node[4] = (int)FaceNodes->GetComponent(face[3], 1);
          node[1] = (int)FaceNodes->GetComponent(face[3], 2);
          node[0] = (int)FaceNodes->GetComponent(face[3], 3);
          }
        }
      else
        {
        BF = face[4];
        if (spinFace[4] > 0)
          {
          node[0] = (int)FaceNodes->GetComponent(face[4], 0);
          node[1] = (int)FaceNodes->GetComponent(face[4], 1);
          node[4] = (int)FaceNodes->GetComponent(face[4], 2);
          node[3] = (int)FaceNodes->GetComponent(face[4], 3);
          }
        else
          {
          node[3] = (int)FaceNodes->GetComponent(face[4], 0);
          node[4] = (int)FaceNodes->GetComponent(face[4], 1);
          node[1] = (int)FaceNodes->GetComponent(face[4], 2);
          node[0] = (int)FaceNodes->GetComponent(face[4], 3);
          }
        }

      int trf[2];
      int index = 0;
      if ( FaceTypes->GetValue(face[0]) == 3)
        {
        trf[index] = face[0];
        index++;
        }
      if ( FaceTypes->GetValue(face[1]) == 3)
        {
        trf[index] = face[1];
        index++;
        }
      if ( FaceTypes->GetValue(face[2]) == 3)
        {
        trf[index] = face[2];
        index++;
        }
      if ( FaceTypes->GetValue(face[3]) == 3)
        {
        trf[index] = face[3];
        index++;
        }
      if ( FaceTypes->GetValue(face[4]) == 3)
        {
        trf[index] = face[4];
        index++;
        }

      tempNode[0] = (int)FaceNodes->GetComponent(trf[0], 0);
      tempNode[1] = (int)FaceNodes->GetComponent(trf[0], 1);
      tempNode[2] = (int)FaceNodes->GetComponent(trf[0], 2);

      tempNode[3] = (int)FaceNodes->GetComponent(trf[1], 0);
      tempNode[4] = (int)FaceNodes->GetComponent(trf[1], 1);
      tempNode[5] = (int)FaceNodes->GetComponent(trf[1], 2);

      if (((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])&&(tempNode[0]!=node[4])&&(tempNode[0]!=node[3]))
        && ((tempNode[1]==node[0])||(tempNode[1]==node[1])) )
        {
        node[2] = tempNode[0];
        }
      else if (((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])&&(tempNode[0]!=node[4])&&(tempNode[0]!=node[3]))
        && ((tempNode[1]==node[4])||(tempNode[1]==node[3])))
        {
        node[5] = tempNode[0];
        }

      if (((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])&&(tempNode[1]!=node[4])&&(tempNode[1]!=node[3]))
        && ((tempNode[0]==node[0])||(tempNode[0]==node[1])))
        {
        node[2] = tempNode[1];
        }
      else if (((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])&&(tempNode[1]!=node[4])&&(tempNode[1]!=node[3]))
        && ((tempNode[0]==node[4])||(tempNode[0]==node[3])))
        {
        node[5] = tempNode[1];
        }

      if (((tempNode[2]!=node[0])&&(tempNode[2]!=node[1])&&(tempNode[2]!=node[4])&&(tempNode[2]!=node[3])) 
        && ((tempNode[1]==node[0])||(tempNode[1]==node[1]))) 
        {
        node[2] = tempNode[2];
        }
      else if (((tempNode[2]!=node[0])&&(tempNode[2]!=node[1])&&(tempNode[2]!=node[4])&&(tempNode[2]!=node[3]))
        && ((tempNode[1]==node[4])||(tempNode[1]==node[3])))
        {
        node[5] = tempNode[2];
        }

      if (((tempNode[3]!=node[0])&&(tempNode[3]!=node[1])&&(tempNode[3]!=node[4])&&(tempNode[3]!=node[3]))
        && ((tempNode[4]==node[0])||(tempNode[4]==node[1])))
        {
        node[2] = tempNode[3];
        }
      else if (((tempNode[3]!=node[0])&&(tempNode[3]!=node[1])&&(tempNode[3]!=node[4])&&(tempNode[3]!=node[3]))
        && ((tempNode[4]==node[4])||(tempNode[4]==node[3])))
        {
        node[5] = tempNode[3];
        }

      if (((tempNode[4]!=node[0])&&(tempNode[4]!=node[1])&&(tempNode[4]!=node[4])&&(tempNode[4]!=node[3])) 
        && ((tempNode[3]==node[0])||(tempNode[3]==node[1])))
        {
        node[2] = tempNode[4];
        }
      else if (((tempNode[4]!=node[0])&&(tempNode[4]!=node[1])&&(tempNode[4]!=node[4])&&(tempNode[4]!=node[3])) 
        && ((tempNode[3]==node[4])||(tempNode[3]==node[3])))
        {
        node[5] = tempNode[4];
        }

      if (((tempNode[5]!=node[0])&&(tempNode[5]!=node[1])&&(tempNode[5]!=node[4])&&(tempNode[5]!=node[3]))
        && ((tempNode[4]==node[0])||(tempNode[4]==node[1])))
        {
        node[2] = tempNode[5];
        }
      else if (((tempNode[5]!=node[0])&&(tempNode[5]!=node[1])&&(tempNode[5]!=node[4])&&(tempNode[5]!=node[3]))
        && ((tempNode[4]==node[4])||(tempNode[4]==node[3])))
        {
        node[5] = tempNode[5];
        }

      AWedge->GetPointIds()->SetId( 0, node[0]);
      AWedge->GetPointIds()->SetId( 1, node[1]);
      AWedge->GetPointIds()->SetId( 2, node[2]);
      AWedge->GetPointIds()->SetId( 3, node[3]);
      AWedge->GetPointIds()->SetId( 4, node[4]);
      AWedge->GetPointIds()->SetId( 5, node[5]);

      if (CellParentFlags->GetValue(i) != 1)
        {
        Mesh->InsertNextCell(AWedge->GetCellType(), AWedge->GetPointIds());
        }
      }
    }
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellParentFlags(void)
{
  // Initialize Array
  for (int i = 1; i <= NumberOfCells; i++)
    {
    CellParentFlags->InsertValue(i,0);
    }

  for (int i = 0; i < NumberOfCellTrees; i++)
    {
    for (int j = CellTreeParentCellId0->GetValue(i); 
      j <= CellTreeParentCellId1->GetValue(i); j++)
      {
      CellParentFlags->InsertValue(j,1);
      }
    }
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellNumberOfFaces(void)
{
  for (int i = 0; i <= NumberOfCells; i++)
    {
    CellNumberOfFaces->InsertValue( i, 0);
    }

  for (int i = 1; i <= NumberOfFaces; i++)
    {
    int c0 = (int)FaceCells->GetComponent( i, 0);
    int c1 = (int)FaceCells->GetComponent( i, 1);
    int nc0 = CellNumberOfFaces->GetValue(c0);
    int nc1 = CellNumberOfFaces->GetValue(c1);

    if ( c0 != 0)
      {
      nc0++;
      CellNumberOfFaces->InsertValue( c0, nc0);
      }

    if ( c1 != 0)
      {
      nc1++;
      CellNumberOfFaces->InsertValue( c1, nc1);
      }
    }
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellFaces(void)
{
  // Make an index array to determine where each cell is in the cell 
  // face array and ...
  // Make a temporary number of faces/cell array to keep track of 
  // where to put the faces within each block.

  int index = 0;
  int *NumberOfFacesInCell;
  NumberOfFacesInCell = new int[NumberOfCells+1];

  for (int i = 1; i <= NumberOfCells; i++)
    {
    CellIndex->InsertValue(i, index);	
    index = index + CellNumberOfFaces->GetValue(i); 
    NumberOfFacesInCell[i] = 0;
    }

  CellIndex->InsertValue(0, 0);
  NumberOfFacesInCell[0] = 0;

  for (int i = 1; i <= NumberOfFaces; i++)
    {
    int c0 = (int)FaceCells->GetComponent(i,0);
    int c1 = (int)FaceCells->GetComponent(i,1);
    int nc0 = NumberOfFacesInCell[c0];
    int nc1 = NumberOfFacesInCell[c1];
    int ic0 = CellIndex->GetValue(c0);
    int ic1 = CellIndex->GetValue(c1);

    if ( c0 != 0)
      {
      CellFaces->InsertValue( ic0+nc0, i);
      nc0++;
      NumberOfFacesInCell[c0] = nc0;
      }

    if ( c1 != 0 )
      {
      CellFaces->InsertValue( ic1 + nc1, i);
      nc1++;
      NumberOfFacesInCell[c1] = nc1;
      }
    }
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::RemoveExtraFaces(void)
{
  int faces[1000000];
  int badKids[1000000];
  int numberOfBadKids = 0;
  int actualFaces[7];

  actualFaces[0] = 0;  // Mixed 
  actualFaces[1] = 3;  // triangular 
  actualFaces[2] = 4;  // tetrahedral
  actualFaces[3] = 4;  // quadrilateral
  actualFaces[4] = 6;  // hexahedral
  actualFaces[5] = 5;  // pyramid
  actualFaces[6] = 5;  // wedge

  // Initialize Clean Cell Array
  for (int i = 0; i <= NumberOfCells; i++)
    {
    for(int j = 0; j < 6; j++)
      {
      CellFacesClean->InsertComponent( i, j, 0);
      }
    }

  for (int i = 1; i <= NumberOfCells; i++)
    {
    numberOfBadKids = 0;
    int cellType = CellTypes->GetValue(i);
    int numberOfFaces = CellNumberOfFaces->GetValue(i);

    if ( numberOfFaces > actualFaces[cellType])
      {
      int ic = CellIndex->GetValue(i);
      for (int j = 0; j < numberOfFaces; j++)
        {
        int face = CellFaces->GetValue(ic+j);
        int parentFlag = FaceParentFlags->GetValue(face);
        int ifChildFlag = InterfaceFaceChildFlags->GetValue(face);
        int ncgFaceChildFlag = NCGFaceChildFlags->GetValue(face);

        if (parentFlag == 1)
          {
          int startKid = 
            FaceTreesKidsIndex->GetValue(FaceTreeParentTable->GetValue(face));
          int endKid =
            FaceTreesKidsIndex->GetValue(FaceTreeParentTable->GetValue(face))
	    + FaceTreesNumberOfKids->
            GetValue(FaceTreeParentTable->GetValue(face));
          for (int k = startKid; k < endKid; k++)
            {
            badKids[numberOfBadKids] = FaceTreesKids->GetValue(k);
            numberOfBadKids++;
            }
          }

        if ( ifChildFlag == 1)
          {
          badKids[numberOfBadKids] = face;
          numberOfBadKids++;
          }

        if (ncgFaceChildFlag == 1)
          {
          badKids[numberOfBadKids] = face;
          numberOfBadKids++;
          }
        }

      if ((numberOfBadKids +actualFaces[cellType]) !=  numberOfFaces)
        {
        cout << " Problem in Face Reduction !!!! " << endl;
        cout << " Cell = " << i << endl;
        cout << " Problem - Number of Faces = " 
          << numberOfFaces << ", Actual Faces " 
          << actualFaces[CellTypes->GetValue(i)] 
          << ", Cell Type = " << CellTypes->GetValue(i) << endl;
        }

      int idx = 0;
      for (int j = 0; j < numberOfFaces; j++)
        {
        int bk = 0;
        int face = CellFaces->GetValue( ic + j);
        for (int m = 0; m < numberOfBadKids; m++)
          {
          if ( badKids[m] == face)
            {
            bk = 1;
            }
          }

          if (bk == 0)
            {
            faces[idx] = face;
            idx++;
            }
        }
      }
    else
      {
      int idx = 0;
      int ic = CellIndex->GetValue(i);
      for (int j = 0; j < numberOfFaces; j++)
        {
        int face = CellFaces->GetValue(ic+j);
        faces[idx] = face;
        idx++;
        }
      }

    for (int j = 0; j < actualFaces[cellType]; j++)
      {
      CellFacesClean->InsertComponent( i, j, faces[j]);
      }
    }
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::ParseDataFile(void)
{
  int bufptr = 0;
  while ( bufptr < DataFileBufferLength)
    {
    if ( DataFileBuffer[bufptr] == '(')
      {
      int ix = GetDataIndex(bufptr);
      bufptr = ExecuteDataTask(ix, bufptr);
      }
    bufptr++;
    }
  return;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::InitializeVariableNames ( void )
{
  VariableNames[1] = "PRESSURE";
  VariableNames[2] = "MOMENTUM";
  VariableNames[3] = "TEMPERATURE";
  VariableNames[4] = "ENTHALPY";
  VariableNames[5] = "TKE";
  VariableNames[6] = "TED";
  VariableNames[7] = "SPECIES";
  VariableNames[8] = "G";
  VariableNames[9] = "WSWIRL";
  VariableNames[10] = "DPMS_MASS";
  VariableNames[11] = "DPMS_MOM";
  VariableNames[12] = "DPMS_ENERGY";
  VariableNames[13] = "DPMS_SPECIES";
  VariableNames[14] = "DVOLUME_DT";
  VariableNames[15] = "BODY_FORCES";
  VariableNames[16] = "FMEAN";
  VariableNames[17] = "FVAR";
  VariableNames[18] = "MASS_FLUX";
  VariableNames[19] = "WALL_SHEAR";
  VariableNames[20] = "BOUNDARY_HEAT_FLUX";
  VariableNames[21] = "BOUNDARY_RAD_HEAT_FLUX";
  VariableNames[22] = "OLD_PRESSURE";
  VariableNames[23] = "POLLUT";
  VariableNames[24] = "DPMS_P1_S";
  VariableNames[25] = "DPMS_P1_AP";
  VariableNames[26] = "WALL_GAS_TEMPERATURE";
  VariableNames[27] = "DPMS_P1_DIFF";
  VariableNames[28] = "DR_SURF";
  VariableNames[29] = "W_M1";
  VariableNames[30] = "W_M2";
  VariableNames[31] = "DPMS_BURNOUT";
  VariableNames[32] = "DPMS_CONCENTRATION";
  VariableNames[33] = "PDF_MW";
  VariableNames[34] = "DPMS_WSWIRL";
  VariableNames[35] = "YPLUS";
  VariableNames[36] = "YPLUS_UTAU";
  VariableNames[37] = "WALL_SHEAR_SWIRL";
  VariableNames[38] = "WALL_T_INNER";
  VariableNames[39] = "POLLUT0";
  VariableNames[40] = "POLLUT1";
  VariableNames[41] = "WALL_G_INNER";
  VariableNames[42] = "PREMIXC";
  VariableNames[43] = "PREMIXC_T";
  VariableNames[44] = "PREMIXC_RATE";
  VariableNames[45] = "POLLUT2";
  VariableNames[46] = "POLLUT3";
  VariableNames[47] = "MASS_FLUX_M1";
  VariableNames[48] = "MASS_FLUX_M2";
  VariableNames[49] = "GRID_FLUX";
  VariableNames[50] = "DO_I";
  VariableNames[51] = "DO_RECON_I";
  VariableNames[52] = "DO_ENERGY_SOURCE";
  VariableNames[53] = "DO_IRRAD";
  VariableNames[54] = "DO_QMINUS";
  VariableNames[55] = "DO_IRRAD_OLD";
  VariableNames[56] = "DO_IWX";
  VariableNames[57] = "DO_IWY";
  VariableNames[58] = "DO_IWZ";
  VariableNames[59] = "MACH";
  VariableNames[60] = "SLIP_U";
  VariableNames[61] = "SLIP_V";
  VariableNames[62] = "SLIP_W";
  VariableNames[63] = "SDR";
  VariableNames[64] = "SDR_M1";
  VariableNames[65] = "SDR_M2";
  VariableNames[66] = "POLLUT4";
  VariableNames[67] = "GRANULAR_TEMPERATURE";
  VariableNames[68] = "GRANULAR_TEMPERATURE_M1";
  VariableNames[69] = "GRANULAR_TEMPERATURE_M2";
  VariableNames[70] = "VFLUX";
  VariableNames[80] = "VFLUX_M1";
  VariableNames[90] = "VFLUX_M2";
  VariableNames[91] = "DO_QNET";
  VariableNames[92] = "DO_QTRANS";
  VariableNames[93] = "DO_QREFL";
  VariableNames[94] = "DO_QABS";
  VariableNames[95] = "POLLUT5";
  VariableNames[96] = "WALL_DIST";
  VariableNames[97] = "SOLAR_SOURCE";
  VariableNames[98] = "SOLAR_QREFL";
  VariableNames[99] = "SOLAR_QABS";
  VariableNames[100] = "SOLAR_QTRANS";
  VariableNames[101] = "DENSITY";
  VariableNames[102] = "MU_LAM";
  VariableNames[103] = "MU_TURB";
  VariableNames[104] = "CP";
  VariableNames[105] = "KTC";
  VariableNames[106] = "VGS_DTRM";
  VariableNames[107] = "VGF_DTRM";
  VariableNames[108] = "RSTRESS";
  VariableNames[109] = "THREAD_RAD_FLUX";
  VariableNames[110] = "SPE_Q";
  VariableNames[111] = "X_VELOCITY";
  VariableNames[112] = "Y_VELOCITY";
  VariableNames[113] = "Z_VELOCITY";
  VariableNames[114] = "WALL_VELOCITY";
  VariableNames[115] = "X_VELOCITY_M1";
  VariableNames[116] = "Y_VELOCITY_M1";
  VariableNames[117] = "Z_VELOCITY_M1";
  VariableNames[118] = "PHASE_MASS";
  VariableNames[119] = "TKE_M1";
  VariableNames[120] = "TED_M1";
  VariableNames[121] = "POLLUT6";
  VariableNames[122] = "X_VELOCITY_M2";
  VariableNames[123] = "Y_VELOCITY_M2";
  VariableNames[124] = "Z_VELOCITY_M2";
  VariableNames[126] = "TKE_M2";
  VariableNames[127] = "TED_M2";
  VariableNames[128] = "RUU";
  VariableNames[129] = "RVV";
  VariableNames[130] = "RWW";
  VariableNames[131] = "RUV";
  VariableNames[132] = "RVW";
  VariableNames[133] = "RUW";
  VariableNames[134] = "DPMS_EROSION";
  VariableNames[135] = "DPMS_ACCRETION";
  VariableNames[136] = "FMEAN2";
  VariableNames[137] = "FVAR2";
  VariableNames[138] = "ENTHALPY_M1";
  VariableNames[139] = "ENTHALPY_M2";
  VariableNames[140] = "FMEAN_M1";
  VariableNames[141] = "FMEAN_M2";
  VariableNames[142] = "FVAR_M1";
  VariableNames[143] = "FVAR_M2";
  VariableNames[144] = "FMEAN2_M1";
  VariableNames[145] = "FMEAN2_M2";
  VariableNames[146] = "FVAR2_M1";
  VariableNames[147] = "FVAR2_M2";
  VariableNames[148] = "PREMIXC_M1";
  VariableNames[149] = "PREMIXC_M2";
  VariableNames[150] = "VOF";
  VariableNames[151] = "VOF_1";
  VariableNames[152] = "VOF_2";
  VariableNames[153] = "VOF_3";
  VariableNames[154] = "VOF_4";
  VariableNames[160] = "VOF_M1";
  VariableNames[161] = "VOF_1_M1";
  VariableNames[162] = "VOF_2_M1";
  VariableNames[163] = "VOF_3_M1";
  VariableNames[164] = "VOF_4_M1";
  VariableNames[170] = "VOF_M2";
  VariableNames[171] = "VOF_1_M2";
  VariableNames[172] = "VOF_2_M2";
  VariableNames[173] = "VOF_3_M2";
  VariableNames[174] = "VOF_4_M2";
  VariableNames[180] = "VOLUME_M2";
  VariableNames[181] = "WALL_GRID_VELOCITY";
  VariableNames[190] = "SV_T_AUX";
  VariableNames[191] = "SV_T_AP_AUX";
  VariableNames[192] = "TOTAL_PRESSURE";
  VariableNames[193] = "TOTAL_TEMPERATURE";
  VariableNames[194] = "NRBC_DC";
  VariableNames[195] = "DP_TMFR";
  VariableNames[200] = "SV_Y_0";
  VariableNames[201] = "SV_Y_1";
  VariableNames[202] = "SV_Y_2";
  VariableNames[203] = "SV_Y_3";
  VariableNames[204] = "SV_Y_4";
  VariableNames[205] = "SV_Y_5";
  VariableNames[206] = "SV_Y_6";
  VariableNames[207] = "SV_Y_7";
  VariableNames[208] = "SV_Y_8";
  VariableNames[209] = "SV_Y_9";
  VariableNames[210] = "SV_Y_10";
  VariableNames[211] = "SV_Y_11";
  VariableNames[212] = "SV_Y_12";
  VariableNames[213] = "SV_Y_13";
  VariableNames[214] = "SV_Y_14";
  VariableNames[215] = "SV_Y_15";
  VariableNames[216] = "SV_Y_16";
  VariableNames[217] = "SV_Y_17";
  VariableNames[218] = "SV_Y_18";
  VariableNames[219] = "SV_Y_19";
  VariableNames[220] = "SV_Y_20";
  VariableNames[221] = "SV_Y_21";
  VariableNames[222] = "SV_Y_22";
  VariableNames[223] = "SV_Y_23";
  VariableNames[224] = "SV_Y_24";
  VariableNames[225] = "SV_Y_25";
  VariableNames[226] = "SV_Y_26";
  VariableNames[227] = "SV_Y_27";
  VariableNames[228] = "SV_Y_28";
  VariableNames[229] = "SV_Y_29";
  VariableNames[230] = "SV_Y_30";
  VariableNames[231] = "SV_Y_31";
  VariableNames[232] = "SV_Y_32";
  VariableNames[233] = "SV_Y_33";
  VariableNames[234] = "SV_Y_34";
  VariableNames[235] = "SV_Y_35";
  VariableNames[236] = "SV_Y_36";
  VariableNames[237] = "SV_Y_37";
  VariableNames[238] = "SV_Y_38";
  VariableNames[239] = "SV_Y_39";
  VariableNames[240] = "SV_Y_40";
  VariableNames[241] = "SV_Y_41";
  VariableNames[242] = "SV_Y_42";
  VariableNames[243] = "SV_Y_43";
  VariableNames[244] = "SV_Y_44";
  VariableNames[245] = "SV_Y_45";
  VariableNames[246] = "SV_Y_46";
  VariableNames[247] = "SV_Y_47";
  VariableNames[248] = "SV_Y_48";
  VariableNames[249] = "SV_Y_49";
  VariableNames[250] = "SV_M1_Y_0";
  VariableNames[251] = "SV_M1_Y_1";
  VariableNames[252] = "SV_M1_Y_2";
  VariableNames[253] = "SV_M1_Y_3";
  VariableNames[254] = "SV_M1_Y_4";
  VariableNames[255] = "SV_M1_Y_5";
  VariableNames[256] = "SV_M1_Y_6";
  VariableNames[257] = "SV_M1_Y_7";
  VariableNames[258] = "SV_M1_Y_8";
  VariableNames[259] = "SV_M1_Y_9";
  VariableNames[260] = "SV_M1_Y_10";
  VariableNames[261] = "SV_M1_Y_11";
  VariableNames[262] = "SV_M1_Y_12";
  VariableNames[263] = "SV_M1_Y_13";
  VariableNames[264] = "SV_M1_Y_14";
  VariableNames[265] = "SV_M1_Y_15";
  VariableNames[266] = "SV_M1_Y_16";
  VariableNames[267] = "SV_M1_Y_17";
  VariableNames[268] = "SV_M1_Y_18";
  VariableNames[269] = "SV_M1_Y_19";
  VariableNames[270] = "SV_M1_Y_20";
  VariableNames[271] = "SV_M1_Y_21";
  VariableNames[272] = "SV_M1_Y_22";
  VariableNames[273] = "SV_M1_Y_23";
  VariableNames[274] = "SV_M1_Y_24";
  VariableNames[275] = "SV_M1_Y_25";
  VariableNames[276] = "SV_M1_Y_26";
  VariableNames[277] = "SV_M1_Y_27";
  VariableNames[278] = "SV_M1_Y_28";
  VariableNames[279] = "SV_M1_Y_29";
  VariableNames[280] = "SV_M1_Y_30";
  VariableNames[281] = "SV_M1_Y_31";
  VariableNames[282] = "SV_M1_Y_32";
  VariableNames[283] = "SV_M1_Y_33";
  VariableNames[284] = "SV_M1_Y_34";
  VariableNames[285] = "SV_M1_Y_35";
  VariableNames[286] = "SV_M1_Y_36";
  VariableNames[287] = "SV_M1_Y_37";
  VariableNames[288] = "SV_M1_Y_38";
  VariableNames[289] = "SV_M1_Y_39";
  VariableNames[290] = "SV_M1_Y_40";
  VariableNames[291] = "SV_M1_Y_41";
  VariableNames[292] = "SV_M1_Y_42";
  VariableNames[293] = "SV_M1_Y_43";
  VariableNames[294] = "SV_M1_Y_44";
  VariableNames[295] = "SV_M1_Y_45";
  VariableNames[296] = "SV_M1_Y_46";
  VariableNames[297] = "SV_M1_Y_47";
  VariableNames[298] = "SV_M1_Y_48";
  VariableNames[299] = "SV_M1_Y_49";
  VariableNames[300] = "SV_M2_Y_0";
  VariableNames[301] = "SV_M2_Y_1";
  VariableNames[302] = "SV_M2_Y_2";
  VariableNames[303] = "SV_M2_Y_3";
  VariableNames[304] = "SV_M2_Y_4";
  VariableNames[305] = "SV_M2_Y_5";
  VariableNames[306] = "SV_M2_Y_6";
  VariableNames[307] = "SV_M2_Y_7";
  VariableNames[308] = "SV_M2_Y_8";
  VariableNames[309] = "SV_M2_Y_9";
  VariableNames[310] = "SV_M2_Y_10";
  VariableNames[311] = "SV_M2_Y_11";
  VariableNames[312] = "SV_M2_Y_12";
  VariableNames[313] = "SV_M2_Y_13";
  VariableNames[314] = "SV_M2_Y_14";
  VariableNames[315] = "SV_M2_Y_15";
  VariableNames[316] = "SV_M2_Y_16";
  VariableNames[317] = "SV_M2_Y_17";
  VariableNames[318] = "SV_M2_Y_18";
  VariableNames[319] = "SV_M2_Y_19";
  VariableNames[320] = "SV_M2_Y_20";
  VariableNames[321] = "SV_M2_Y_21";
  VariableNames[322] = "SV_M2_Y_22";
  VariableNames[323] = "SV_M2_Y_23";
  VariableNames[324] = "SV_M2_Y_24";
  VariableNames[325] = "SV_M2_Y_25";
  VariableNames[326] = "SV_M2_Y_26";
  VariableNames[327] = "SV_M2_Y_27";
  VariableNames[328] = "SV_M2_Y_28";
  VariableNames[329] = "SV_M2_Y_29";
  VariableNames[330] = "SV_M2_Y_30";
  VariableNames[331] = "SV_M2_Y_31";
  VariableNames[332] = "SV_M2_Y_32";
  VariableNames[333] = "SV_M2_Y_33";
  VariableNames[334] = "SV_M2_Y_34";
  VariableNames[335] = "SV_M2_Y_35";
  VariableNames[336] = "SV_M2_Y_36";
  VariableNames[337] = "SV_M2_Y_37";
  VariableNames[338] = "SV_M2_Y_38";
  VariableNames[339] = "SV_M2_Y_39";
  VariableNames[340] = "SV_M2_Y_40";
  VariableNames[341] = "SV_M2_Y_41";
  VariableNames[342] = "SV_M2_Y_42";
  VariableNames[343] = "SV_M2_Y_43";
  VariableNames[344] = "SV_M2_Y_44";
  VariableNames[345] = "SV_M2_Y_45";
  VariableNames[346] = "SV_M2_Y_46";
  VariableNames[347] = "SV_M2_Y_47";
  VariableNames[348] = "SV_M2_Y_48";
  VariableNames[349] = "SV_M2_Y_49";
  VariableNames[350] = "DR_SURF_0";
  VariableNames[351] = "DR_SURF_1";
  VariableNames[352] = "DR_SURF_2";
  VariableNames[353] = "DR_SURF_3";
  VariableNames[354] = "DR_SURF_4";
  VariableNames[355] = "DR_SURF_5";
  VariableNames[356] = "DR_SURF_6";
  VariableNames[357] = "DR_SURF_7";
  VariableNames[358] = "DR_SURF_8";
  VariableNames[359] = "DR_SURF_9";
  VariableNames[360] = "DR_SURF_10";
  VariableNames[361] = "DR_SURF_11";
  VariableNames[362] = "DR_SURF_12";
  VariableNames[363] = "DR_SURF_13";
  VariableNames[364] = "DR_SURF_14";
  VariableNames[365] = "DR_SURF_15";
  VariableNames[366] = "DR_SURF_16";
  VariableNames[367] = "DR_SURF_17";
  VariableNames[368] = "DR_SURF_18";
  VariableNames[369] = "DR_SURF_19";
  VariableNames[370] = "DR_SURF_20";
  VariableNames[371] = "DR_SURF_21";
  VariableNames[372] = "DR_SURF_22";
  VariableNames[373] = "DR_SURF_23";
  VariableNames[374] = "DR_SURF_24";
  VariableNames[375] = "DR_SURF_25";
  VariableNames[376] = "DR_SURF_26";
  VariableNames[377] = "DR_SURF_27";
  VariableNames[378] = "DR_SURF_28";
  VariableNames[379] = "DR_SURF_29";
  VariableNames[380] = "DR_SURF_30";
  VariableNames[381] = "DR_SURF_31";
  VariableNames[382] = "DR_SURF_32";
  VariableNames[383] = "DR_SURF_33";
  VariableNames[384] = "DR_SURF_34";
  VariableNames[385] = "DR_SURF_35";
  VariableNames[386] = "DR_SURF_36";
  VariableNames[387] = "DR_SURF_37";
  VariableNames[388] = "DR_SURF_38";
  VariableNames[389] = "DR_SURF_39";
  VariableNames[390] = "DR_SURF_40";
  VariableNames[391] = "DR_SURF_41";
  VariableNames[392] = "DR_SURF_42";
  VariableNames[393] = "DR_SURF_43";
  VariableNames[394] = "DR_SURF_44";
  VariableNames[395] = "DR_SURF_45";
  VariableNames[396] = "DR_SURF_46";
  VariableNames[397] = "DR_SURF_47";
  VariableNames[398] = "DR_SURF_48";
  VariableNames[399] = "DR_SURF_49";
  VariableNames[400] = "PRESSURE_MEAN";
  VariableNames[401] = "PRESSURE_RMS";
  VariableNames[402] = "X_VELOCITY_MEAN";
  VariableNames[403] = "X_VELOCITY_RMS";
  VariableNames[404] = "Y_VELOCITY_MEAN";
  VariableNames[405] = "Y_VELOCITY_RMS";
  VariableNames[406] = "Z_VELOCITY_MEAN";
  VariableNames[407] = "Z_VELOCITY_RMS";
  VariableNames[408] = "TEMPERATURE_MEAN";
  VariableNames[409] = "TEMPERATURE_RMS";
  VariableNames[410] = "VOF_MEAN";
  VariableNames[411] = "VOF_RMS";
  VariableNames[412] = "PRESSURE_M1";
  VariableNames[413] = "PRESSURE_M2";
  VariableNames[414] = "GRANULAR_TEMPERATURE_MEAN";
  VariableNames[415] = "GRANULAR_TEMPERATURE_RMS";
  VariableNames[450] = "DPMS_Y_0";
  VariableNames[451] = "DPMS_Y_1";
  VariableNames[452] = "DPMS_Y_2";
  VariableNames[453] = "DPMS_Y_3";
  VariableNames[454] = "DPMS_Y_4";
  VariableNames[455] = "DPMS_Y_5";
  VariableNames[456] = "DPMS_Y_6";
  VariableNames[457] = "DPMS_Y_7";
  VariableNames[458] = "DPMS_Y_8";
  VariableNames[459] = "DPMS_Y_9";
  VariableNames[460] = "DPMS_Y_10";
  VariableNames[461] = "DPMS_Y_11";
  VariableNames[462] = "DPMS_Y_12";
  VariableNames[463] = "DPMS_Y_13";
  VariableNames[464] = "DPMS_Y_14";
  VariableNames[465] = "DPMS_Y_15";
  VariableNames[466] = "DPMS_Y_16";
  VariableNames[467] = "DPMS_Y_17";
  VariableNames[468] = "DPMS_Y_18";
  VariableNames[469] = "DPMS_Y_19";
  VariableNames[470] = "DPMS_Y_20";
  VariableNames[471] = "DPMS_Y_21";
  VariableNames[472] = "DPMS_Y_22";
  VariableNames[473] = "DPMS_Y_23";
  VariableNames[474] = "DPMS_Y_24";
  VariableNames[475] = "DPMS_Y_25";
  VariableNames[476] = "DPMS_Y_26";
  VariableNames[477] = "DPMS_Y_27";
  VariableNames[478] = "DPMS_Y_28";
  VariableNames[479] = "DPMS_Y_29";
  VariableNames[480] = "DPMS_Y_30";
  VariableNames[481] = "DPMS_Y_31";
  VariableNames[482] = "DPMS_Y_32";
  VariableNames[483] = "DPMS_Y_33";
  VariableNames[484] = "DPMS_Y_34";
  VariableNames[485] = "DPMS_Y_35";
  VariableNames[486] = "DPMS_Y_36";
  VariableNames[487] = "DPMS_Y_37";
  VariableNames[488] = "DPMS_Y_38";
  VariableNames[489] = "DPMS_Y_39";
  VariableNames[490] = "DPMS_Y_40";
  VariableNames[491] = "DPMS_Y_41";
  VariableNames[492] = "DPMS_Y_42";
  VariableNames[493] = "DPMS_Y_43";
  VariableNames[494] = "DPMS_Y_44";
  VariableNames[495] = "DPMS_Y_45";
  VariableNames[496] = "DPMS_Y_46";
  VariableNames[497] = "DPMS_Y_47";
  VariableNames[498] = "DPMS_Y_48";
  VariableNames[499] = "DPMS_Y_49";
  VariableNames[500] = "NUT";
  VariableNames[501] = "NUT_M1";
  VariableNames[502] = "NUT_M2";
  VariableNames[503] = "RUU_M1";
  VariableNames[504] = "RVV_M1";
  VariableNames[505] = "RWW_M1";
  VariableNames[506] = "RUV_M1";
  VariableNames[507] = "RVW_M1";
  VariableNames[508] = "RUW_M1";
  VariableNames[509] = "RUU_M2";
  VariableNames[510] = "RVV_M2";
  VariableNames[511] = "RWW_M2";
  VariableNames[512] = "RUV_M2";
  VariableNames[513] = "RVW_M2";
  VariableNames[514] = "RUW_M2";
  VariableNames[515] = "ENERGY_M1";
  VariableNames[516] = "ENERGY_M2";
  VariableNames[517] = "DENSITY_M1";
  VariableNames[518] = "DENSITY_M2";
  VariableNames[519] = "DPMS_PDF_1";
  VariableNames[520] = "DPMS_PDF_2";
  VariableNames[521] = "V2";
  VariableNames[522] = "V2_M1";
  VariableNames[523] = "V2_M2";
  VariableNames[524] = "FEL";
  VariableNames[525] = "FEL_M1";
  VariableNames[526] = "FEL_M2";
  VariableNames[530] = "SHELL_CELL_T";
  VariableNames[531] = "SHELL_FACE_T";
  VariableNames[540] = "DPMS_TKE";
  VariableNames[541] = "DPMS_D";
  VariableNames[542] = "DPMS_O";
  VariableNames[543] = "DPMS_TKE_RUU";
  VariableNames[544] = "DPMS_TKE_RVV";
  VariableNames[545] = "DPMS_TKE_RWW";
  VariableNames[546] = "DPMS_TKE_RUV";
  VariableNames[547] = "DPMS_TKE_RVW";
  VariableNames[548] = "DPMS_TKE_RUW";
  VariableNames[549] = "DPMS_DS_MASS";
  VariableNames[550] = "DPMS_DS_ENERGY";
  VariableNames[551] = "DPMS_DS_TKE";
  VariableNames[552] = "DPMS_DS_D";
  VariableNames[553] = "DPMS_DS_O";
  VariableNames[554] = "DPMS_DS_TKE_RUU";
  VariableNames[555] = "DPMS_DS_TKE_RVV";
  VariableNames[556] = "DPMS_DS_TKE_RWW";
  VariableNames[557] = "DPMS_DS_TKE_RUV";
  VariableNames[558] = "DPMS_DS_TKE_RVW";
  VariableNames[559] = "DPMS_DS_TKE_RUW";
  VariableNames[560] = "DPMS_DS_PDF_1";
  VariableNames[561] = "DPMS_DS_PDF_2";
  VariableNames[562] = "DPMS_DS_EMISS";
  VariableNames[563] = "DPMS_DS_ABS";
  VariableNames[564] = "DPMS_DS_SCAT";
  VariableNames[565] = "DPMS_DS_BURNOUT";
  VariableNames[566] = "DPMS_DS_MOM";
  VariableNames[567] = "DPMS_DS_WSWIRL";
  VariableNames[600] = "DELH";
  VariableNames[601] = "DPMS_MOM_AP";
  VariableNames[602] = "DPMS_WSWIRL_AP";
  VariableNames[603] = "X_PULL";
  VariableNames[604] = "Y_PULL";
  VariableNames[605] = "Z_PULL";
  VariableNames[606] = "LIQF";
  VariableNames[610] = "PDFT_QBAR";
  VariableNames[611] = "PDFT_PHI";
  VariableNames[612] = "PDFT_Q_TA";
  VariableNames[613] = "PDFT_SVOL_TA";
  VariableNames[614] = "PDFT_MASS_TA";
  VariableNames[615] = "PDFT_T4_TA";
  VariableNames[630] = "SCAD_LES";
  VariableNames[645] = "CREV_MASS";
  VariableNames[646] = "CREV_ENRG";
  VariableNames[647] = "CREV_MOM";
  VariableNames[650] = "XF_ACOUSTICS_MODEL";
  VariableNames[651] = "XF_RF_AC_RECEIVERS_DATA";
  VariableNames[652] = "SV_DPDT_RMS";
  VariableNames[653] = "SV_PRESSURE_M1";
  VariableNames[654] = "XF_RF_AC_PERIODIC_INDEX";
  VariableNames[655] = "XF_RF_AC_PERIODIC_PS";
  VariableNames[656] = "XF_RF_AC_F_NORMAL";
  VariableNames[657] = "XF_RF_AC_F_CENTROID";
  VariableNames[660] = "IGNITE";
  VariableNames[661] = "IGNITE_M1";
  VariableNames[662] = "IGNITE_M2";
  VariableNames[663] = "IGNITE_RATE";
  VariableNames[700] = "UDS_0";
  VariableNames[701] = "UDS_1";
  VariableNames[702] = "UDS_2";
  VariableNames[703] = "UDS_3";
  VariableNames[704] = "UDS_4";
  VariableNames[705] = "UDS_5";
  VariableNames[706] = "UDS_6";
  VariableNames[707] = "UDS_7";
  VariableNames[708] = "UDS_8";
  VariableNames[709] = "UDS_9";
  VariableNames[710] = "UDS_10";
  VariableNames[711] = "UDS_11";
  VariableNames[712] = "UDS_12";
  VariableNames[713] = "UDS_13";
  VariableNames[714] = "UDS_14";
  VariableNames[715] = "UDS_15";
  VariableNames[716] = "UDS_16";
  VariableNames[717] = "UDS_17";
  VariableNames[718] = "UDS_18";
  VariableNames[719] = "UDS_19";
  VariableNames[720] = "UDS_20";
  VariableNames[721] = "UDS_21";
  VariableNames[722] = "UDS_22";
  VariableNames[723] = "UDS_23";
  VariableNames[724] = "UDS_24";
  VariableNames[725] = "UDS_25";
  VariableNames[726] = "UDS_26";
  VariableNames[727] = "UDS_27";
  VariableNames[728] = "UDS_28";
  VariableNames[729] = "UDS_29";
  VariableNames[730] = "UDS_30";
  VariableNames[731] = "UDS_31";
  VariableNames[732] = "UDS_32";
  VariableNames[733] = "UDS_33";
  VariableNames[734] = "UDS_34";
  VariableNames[735] = "UDS_35";
  VariableNames[736] = "UDS_36";
  VariableNames[737] = "UDS_37";
  VariableNames[738] = "UDS_38";
  VariableNames[739] = "UDS_39";
  VariableNames[740] = "UDS_40";
  VariableNames[741] = "UDS_41";
  VariableNames[742] = "UDS_42";
  VariableNames[743] = "UDS_43";
  VariableNames[744] = "UDS_44";
  VariableNames[745] = "UDS_45";
  VariableNames[746] = "UDS_46";
  VariableNames[747] = "UDS_47";
  VariableNames[748] = "UDS_48";
  VariableNames[749] = "UDS_49";
  VariableNames[750] = "UDS_M1_0";
  VariableNames[751] = "UDS_M1_1";
  VariableNames[752] = "UDS_M1_2";
  VariableNames[753] = "UDS_M1_3";
  VariableNames[754] = "UDS_M1_4";
  VariableNames[755] = "UDS_M1_5";
  VariableNames[756] = "UDS_M1_6";
  VariableNames[757] = "UDS_M1_7";
  VariableNames[758] = "UDS_M1_8";
  VariableNames[759] = "UDS_M1_9";
  VariableNames[760] = "UDS_M1_10";
  VariableNames[761] = "UDS_M1_11";
  VariableNames[762] = "UDS_M1_12";
  VariableNames[763] = "UDS_M1_13";
  VariableNames[764] = "UDS_M1_14";
  VariableNames[765] = "UDS_M1_15";
  VariableNames[766] = "UDS_M1_16";
  VariableNames[767] = "UDS_M1_17";
  VariableNames[768] = "UDS_M1_18";
  VariableNames[769] = "UDS_M1_19";
  VariableNames[770] = "UDS_M1_20";
  VariableNames[771] = "UDS_M1_21";
  VariableNames[772] = "UDS_M1_22";
  VariableNames[773] = "UDS_M1_23";
  VariableNames[774] = "UDS_M1_24";
  VariableNames[775] = "UDS_M1_25";
  VariableNames[776] = "UDS_M1_26";
  VariableNames[777] = "UDS_M1_27";
  VariableNames[778] = "UDS_M1_28";
  VariableNames[779] = "UDS_M1_29";
  VariableNames[780] = "UDS_M1_30";
  VariableNames[781] = "UDS_M1_31";
  VariableNames[782] = "UDS_M1_32";
  VariableNames[783] = "UDS_M1_33";
  VariableNames[784] = "UDS_M1_34";
  VariableNames[785] = "UDS_M1_35";
  VariableNames[786] = "UDS_M1_36";
  VariableNames[787] = "UDS_M1_37";
  VariableNames[788] = "UDS_M1_38";
  VariableNames[789] = "UDS_M1_39";
  VariableNames[790] = "UDS_M1_40";
  VariableNames[791] = "UDS_M1_41";
  VariableNames[792] = "UDS_M1_42";
  VariableNames[793] = "UDS_M1_43";
  VariableNames[794] = "UDS_M1_44";
  VariableNames[795] = "UDS_M1_45";
  VariableNames[796] = "UDS_M1_46";
  VariableNames[797] = "UDS_M1_47";
  VariableNames[798] = "UDS_M1_48";
  VariableNames[799] = "UDS_M1_49";
  VariableNames[800] = "UDS_M2_0";
  VariableNames[801] = "UDS_M2_1";
  VariableNames[802] = "UDS_M2_2";
  VariableNames[803] = "UDS_M2_3";
  VariableNames[804] = "UDS_M2_4";
  VariableNames[805] = "UDS_M2_5";
  VariableNames[806] = "UDS_M2_6";
  VariableNames[807] = "UDS_M2_7";
  VariableNames[808] = "UDS_M2_8";
  VariableNames[809] = "UDS_M2_9";
  VariableNames[810] = "UDS_M2_10";
  VariableNames[811] = "UDS_M2_11";
  VariableNames[812] = "UDS_M2_12";
  VariableNames[813] = "UDS_M2_13";
  VariableNames[814] = "UDS_M2_14";
  VariableNames[815] = "UDS_M2_15";
  VariableNames[816] = "UDS_M2_16";
  VariableNames[817] = "UDS_M2_17";
  VariableNames[818] = "UDS_M2_18";
  VariableNames[819] = "UDS_M2_19";
  VariableNames[820] = "UDS_M2_20";
  VariableNames[821] = "UDS_M2_21";
  VariableNames[822] = "UDS_M2_22";
  VariableNames[823] = "UDS_M2_23";
  VariableNames[824] = "UDS_M2_24";
  VariableNames[825] = "UDS_M2_25";
  VariableNames[826] = "UDS_M2_26";
  VariableNames[827] = "UDS_M2_27";
  VariableNames[828] = "UDS_M2_28";
  VariableNames[829] = "UDS_M2_29";
  VariableNames[830] = "UDS_M2_30";
  VariableNames[831] = "UDS_M2_31";
  VariableNames[832] = "UDS_M2_32";
  VariableNames[833] = "UDS_M2_33";
  VariableNames[834] = "UDS_M2_34";
  VariableNames[835] = "UDS_M2_35";
  VariableNames[836] = "UDS_M2_36";
  VariableNames[837] = "UDS_M2_37";
  VariableNames[838] = "UDS_M2_38";
  VariableNames[839] = "UDS_M2_39";
  VariableNames[840] = "UDS_M2_40";
  VariableNames[841] = "UDS_M2_41";
  VariableNames[842] = "UDS_M2_42";
  VariableNames[843] = "UDS_M2_43";
  VariableNames[844] = "UDS_M2_44";
  VariableNames[845] = "UDS_M2_45";
  VariableNames[846] = "UDS_M2_46";
  VariableNames[847] = "UDS_M2_47";
  VariableNames[848] = "UDS_M2_48";
  VariableNames[849] = "UDS_M2_49";
  VariableNames[850] = "DPMS_DS_Y_0";
  VariableNames[851] = "DPMS_DS_Y_1";
  VariableNames[852] = "DPMS_DS_Y_2";
  VariableNames[853] = "DPMS_DS_Y_3";
  VariableNames[854] = "DPMS_DS_Y_4";
  VariableNames[855] = "DPMS_DS_Y_5";
  VariableNames[856] = "DPMS_DS_Y_6";
  VariableNames[857] = "DPMS_DS_Y_7";
  VariableNames[858] = "DPMS_DS_Y_8";
  VariableNames[859] = "DPMS_DS_Y_9";
  VariableNames[860] = "DPMS_DS_Y_10";
  VariableNames[861] = "DPMS_DS_Y_11";
  VariableNames[862] = "DPMS_DS_Y_12";
  VariableNames[863] = "DPMS_DS_Y_13";
  VariableNames[864] = "DPMS_DS_Y_14";
  VariableNames[865] = "DPMS_DS_Y_15";
  VariableNames[866] = "DPMS_DS_Y_16";
  VariableNames[867] = "DPMS_DS_Y_17";
  VariableNames[868] = "DPMS_DS_Y_18";
  VariableNames[869] = "DPMS_DS_Y_19";
  VariableNames[870] = "DPMS_DS_Y_20";
  VariableNames[871] = "DPMS_DS_Y_21";
  VariableNames[872] = "DPMS_DS_Y_22";
  VariableNames[873] = "DPMS_DS_Y_23";
  VariableNames[874] = "DPMS_DS_Y_24";
  VariableNames[875] = "DPMS_DS_Y_25";
  VariableNames[876] = "DPMS_DS_Y_26";
  VariableNames[877] = "DPMS_DS_Y_27";
  VariableNames[878] = "DPMS_DS_Y_28";
  VariableNames[879] = "DPMS_DS_Y_29";
  VariableNames[880] = "DPMS_DS_Y_30";
  VariableNames[881] = "DPMS_DS_Y_31";
  VariableNames[882] = "DPMS_DS_Y_32";
  VariableNames[883] = "DPMS_DS_Y_33";
  VariableNames[884] = "DPMS_DS_Y_34";
  VariableNames[885] = "DPMS_DS_Y_35";
  VariableNames[886] = "DPMS_DS_Y_36";
  VariableNames[887] = "DPMS_DS_Y_37";
  VariableNames[888] = "DPMS_DS_Y_38";
  VariableNames[889] = "DPMS_DS_Y_39";
  VariableNames[890] = "DPMS_DS_Y_40";
  VariableNames[891] = "DPMS_DS_Y_41";
  VariableNames[892] = "DPMS_DS_Y_42";
  VariableNames[893] = "DPMS_DS_Y_43";
  VariableNames[894] = "DPMS_DS_Y_44";
  VariableNames[895] = "DPMS_DS_Y_45";
  VariableNames[896] = "DPMS_DS_Y_46";
  VariableNames[897] = "DPMS_DS_Y_47";
  VariableNames[898] = "DPMS_DS_Y_48";
  VariableNames[899] = "DPMS_DS_Y_49";
  VariableNames[910] = "GRANULAR_PRESSURE";
  VariableNames[911] = "DPMS_DS_P1_S";
  VariableNames[912] = "DPMS_DS_P1_AP";
  VariableNames[913] = "DPMS_DS_P1_DIFF";
  VariableNames[970] = "UDM_I";
  VariableNames[1301] = "WSB";
  VariableNames[1302] = "WSN";
  VariableNames[1303] = "WSR";
  VariableNames[1304] = "WSB_M1";
  VariableNames[1305] = "WSB_M2";
  VariableNames[1306] = "WSN_M1";
  VariableNames[1307] = "WSN_M2";
  VariableNames[1308] = "WSR_M1";
  VariableNames[1309] = "WSR_M2";
  VariableNames[1310] = "MASGEN";
  VariableNames[1311] = "NUCRAT";
  VariableNames[1330] = "TEMPERATURE_M1";
  VariableNames[1331] = "TEMPERATURE_M2";
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCaseIndex(int ix)
{
  char b[5];

  if ( CaseFileBuffer[ix+2] == ' ')
    {
    b[0] = CaseFileBuffer[ix+1];
    b[1] = 0;
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( CaseFileBuffer[ix+3] == ' ')
    {
    b[0] = CaseFileBuffer[ix+1];
    b[1] = CaseFileBuffer[ix+2];
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( CaseFileBuffer[ix+4] == ' ')
    {
    b[0] = CaseFileBuffer[ix+1];
    b[1] = CaseFileBuffer[ix+2];
    b[2] = CaseFileBuffer[ix+3];
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if (CaseFileBuffer[ix+5] == ' ')
    {
    b[0] = CaseFileBuffer[ix+1];
    b[1] = CaseFileBuffer[ix+2];
    b[2] = CaseFileBuffer[ix+3];
    b[3] = CaseFileBuffer[ix+4];
    b[4] = 0;
    return atoi(b);
    }
  else
    {
    b[0] = CaseFileBuffer[ix+1];
    b[1] = CaseFileBuffer[ix+2];
    b[2] = CaseFileBuffer[ix+3];
    b[3] = CaseFileBuffer[ix+4];
    b[4] = 0;
    return -1; 
    }
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::ExecuteCaseTask(int task, int file_index)
{
  int new_index = 0;

  switch ( task )
    {
    //
    //  ASCII Area
    //
    case 0:
      new_index = GetComment(file_index);
      break;
    case 1:
      new_index = GetHeader(file_index);
      break;
    case 2:
      new_index = GetGridDimension(file_index);
      break;
    case 4:
      new_index = GetMachineConfiguration(file_index);
      break;
    case 10:
      new_index = GetNodesASCII(file_index);
      break;
    case 12:
      new_index = GetCellsASCII(file_index);
      break;
    case 13:
      new_index = GetFacesASCII(file_index);
      break;
    case 18:
      new_index = GetPeriodicShadowFacesASCII(file_index);
      break;
    case 33:
      new_index = GetGridSizeASCII(file_index);
      break;
    case 37:
      new_index = GetVariablesASCII(file_index);
      break;
    case 38:
      new_index = GetCortexVariablesASCII(file_index);
      break;
    case 39:
      new_index = GetZoneSectionsASCII(file_index);
      break;
    case 40:
      new_index = GetPartitionASCII(file_index);
      break;
    case 41:
      new_index = GetNodeFlagsASCII(file_index);
      break;
    case 45:
      new_index = GetZoneSectionsASCII(file_index);
      break;
    case 54:
      new_index = Command54(file_index);
      break;
    case 58:
      new_index = GetCellTreeASCII(file_index);
      break;
    case 59:
      new_index = GetFaceTreeASCII(file_index);
      break;
    case 61:
      new_index = GetFaceParentsASCII(file_index);
      break;
    case 62:
      new_index = GetNCG1InformationASCII(file_index);
      break;
    case 63:
      new_index = GetNCG2InformationASCII(file_index);
      break;
    case 64:
      new_index = GetDomainVariablesASCII(file_index);
      break;

  //
  // Single Precision
  //

    case 2010:
      new_index = GetNodesSinglePrecision(file_index);
      break;
    case 2012:
      new_index = GetCellsSinglePrecision(file_index);
      break;
    case 2013:
      new_index = GetFacesSinglePrecision(file_index);
      break;
    case 2018:
      new_index = GetPeriodicShadowFacesSinglePrecision(file_index);
      break;
    case 2033:
      new_index = GetGridSizeSinglePrecision(file_index);
      break;
    case 2037:
      new_index = GetVariablesSinglePrecision(file_index);
      break;
    case 2038:
      new_index = GetCortexVariablesSinglePrecision(file_index);
      break;
    case 2039:
      new_index = GetZoneSectionsSinglePrecision(file_index);
      break;
    case 2040:
      new_index = GetPartitionSinglePrecision(file_index);
      break;
    case 2041:
      new_index = GetNodeFlagsSinglePrecision(file_index);
      break;
    case 2045:
      new_index = GetZoneSectionsSinglePrecision(file_index);
      break;
    case 2058:
      new_index = GetCellTreeSinglePrecision(file_index);
      break;
    case 2059:
      new_index = GetFaceTreeSinglePrecision(file_index);
      break;
    case 2061:
      new_index = GetFaceParentsSinglePrecision(file_index);
      break;
    case 2062:
      new_index = GetNCG1InformationSinglePrecision(file_index);
      break;
    case 2063:
      new_index = GetNCG2InformationSinglePrecision(file_index);
      break;
    case 2064:
      new_index = GetDomainVariablesSinglePrecision(file_index);
      break;

//
// Double Precision
//

    case 3010:
      new_index = GetNodesDoublePrecision(file_index);
      break;
    case 3012:
      new_index = GetCellsDoublePrecision(file_index);
      break;
    case 3013:
      new_index = GetFacesDoublePrecision(file_index);
      break;
    case 3018:
      new_index = GetPeriodicShadowFacesDoublePrecision(file_index);
      break;
    case 3033:
      new_index = GetGridSizeDoublePrecision(file_index);
      break;
    case 3037:
      new_index = GetVariablesDoublePrecision(file_index);
      break;
    case 3038:
      new_index = GetCortexVariablesDoublePrecision(file_index);
      break;
    case 3039:
      new_index = GetZoneSectionsDoublePrecision(file_index);
      break;
    case 3040:
      new_index = GetPartitionDoublePrecision(file_index);
      break;
    case 3041:
      new_index = GetNodeFlagsDoublePrecision(file_index);
      break;
    case 3045:
      new_index = GetZoneSectionsDoublePrecision(file_index);
      break;
    case 3058:
      new_index = GetCellTreeDoublePrecision(file_index);
      break;
    case 3059:
      new_index = GetFaceTreeDoublePrecision(file_index);
      break;
    case 3061:
      new_index = GetFaceParentsDoublePrecision(file_index);
      break;
    case 3062:
      new_index = GetNCG1InformationDoublePrecision(file_index);
      break;
    case 3063:
      new_index = GetNCG2InformationDoublePrecision(file_index);
      break;
    case 3064:
      new_index = GetDomainVariablesDoublePrecision(file_index);
      break;
    default:
      cout << " Unknown Index " << task << endl;
      break;
    }

  return new_index;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataIndex(int ix)
{
  char b[5];

  if ( DataFileBuffer[ix+2] == ' ')
    {
    b[0] = DataFileBuffer[ix+1];
    b[1] = 0;
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( DataFileBuffer[ix+3] == ' ')
    {
    b[0] = DataFileBuffer[ix+1];
    b[1] = DataFileBuffer[ix+2];
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( DataFileBuffer[ix+4] == ' ')
    {
    b[0] = DataFileBuffer[ix+1];
    b[1] = DataFileBuffer[ix+2];
    b[2] = DataFileBuffer[ix+3];
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( DataFileBuffer[ix+5] == ' ')
    {
    b[0] = DataFileBuffer[ix+1];
    b[1] = DataFileBuffer[ix+2];
    b[2] = DataFileBuffer[ix+3];
    b[3] = DataFileBuffer[ix+4];
    b[4] = 0;
    return atoi(b);
    }
  else
    {
    b[0] = DataFileBuffer[ix+1];
    b[1] = DataFileBuffer[ix+2];
    b[2] = DataFileBuffer[ix+3];
    b[3] = DataFileBuffer[ix+4];
    b[4] = 0;
    return -1; 
    }
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::ExecuteDataTask(int task, int file_index)
{
  int new_index;

  switch ( task )
    {
    case 0:
      new_index = GetDataComment(file_index);
      break;
    case 1:
      new_index = GetDataHeader(file_index);
      break;
    case 2:
      new_index = GetDataGridDimension(file_index);
      break;
    case 4:
      new_index = GetDataMachineConfiguration(file_index);
      break;
    case 33:
      new_index = GetDataGridSizeASCII(file_index);
      break;
    case 37:
      new_index = GetDataVariablesASCII(file_index);
      break;
    case 38:
      new_index = GetDataVariablesASCII(file_index);
      break;
    case 50:
      new_index = GetDataUnknownASCII(file_index);
      break;
    case 64:
      new_index = GetDataVariablesASCII(file_index);
      break;
    case 300:
      new_index = GetDataASCII(file_index);
      break;
    case 301:
      new_index = GetUnknownASCII301(file_index);
      break;
    case 302:
      new_index = GetUnknownASCII302(file_index);
      break;
    case 303:
      new_index = GetUnknownASCII303(file_index);
      break;
    case 313:
      new_index = GetUnknownASCII313(file_index);
      break;
    case 2300:
      new_index = GetDataSinglePrecision(file_index);
      break;
    case 2301:
      new_index = GetUnknownSinglePrecision2301(file_index);
      break;
    case 2302:
      new_index = GetUnknownSinglePrecision2302(file_index);
      break;
    case 2303:
      new_index = GetUnknownSinglePrecision2302(file_index);
      break;
    case 2313:
      new_index = GetUnknownSinglePrecision2313(file_index);
      break;
    case 3300:
      new_index = GetDataDoublePrecision(file_index);
      break;
    case 3301:
      new_index = GetUnknownDoublePrecision3301(file_index);
      break;
    case 3302:
      new_index = GetUnknownDoublePrecision3302(file_index);
      break;
    case 3313:
      new_index = GetUnknownDoublePrecision3313(file_index);
      break;
    default:
      cout << " Unknown Index " << task << endl;
      exit(1);
      break;
    }
  return new_index;
}


//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetComment(int ix)
{
  return GoToNextRightParen(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetHeader(int ix)
{
  return GoToNextRightParen(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetMachineConfiguration(int ix)
{
  char buf[120];
  int j = ix+1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );
  j = GoToNextRightParen(j)+1;

  int a, b, c, d, e, f, g, h, m, n, o;
  sscanf( buf, " %d %d %d %d %d %d %d %d %d %d %d", 
    &a, &b, &c, &d, &e, &f, &g, &h, &m, &n, &o );

  if ( a == 60 )
    {
    LittleEndianFlag = 1;
    }
  else
    {
    LittleEndianFlag = 0;
    }

  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetVariablesASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCortexVariablesASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDomainVariablesASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsASCII(int ix)
{
  char buf[120];
  int j = ix+1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );
  j = GoToNextRightParen(j)+1;

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if ( zi != 0)
    {
    CellZones->InsertValue(NumberOfCellZones, zi);
    NumberOfCellZones++;
    }

  if ( zi == 0) 
    {
    NumberOfCells = li;
    }
  else
    {
    if ( et == 0)
      {
      GetMixedCellTypes( j, fi, li); 
      }
    else
      {
      for (int i = fi; i <= li; i++)
        {
        CellTypes->InsertValue(i, et);
        }
      }
    }
  return GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if (zi == 0)
    {
    NumberOfFaces = li;
    }
  else
    {
    j = GoToNextLeftParen(j)+1;
    j = GoToNextEOL(j) +1;
    int n0, n1, n2, n3;
    int c0, c1;
    int type;
    for (int k = fi; k <= li; k++)
      {
      GetStringToNextRightParenOrEOL( j, buf );
      if ( et == 0 )
        {
        if ( buf[0] == 2)
          {
          sscanf( buf, " %x %x %x %x %x ", &type, &n0 , &n1, &c0, &c1 );
          FaceTypes->InsertValue(k,type);
          FaceNodes->InsertComponent(k,0,n0);
          FaceNodes->InsertComponent(k,1,n1);
          FaceNodes->InsertComponent(k,2,0);
          FaceNodes->InsertComponent(k,3,0);
          FaceCells->InsertComponent(k,0,c0);
          FaceCells->InsertComponent(k,1,c1);
          }
        else if ( buf[1] == 3)
          {
          sscanf( buf, " %x %x %x %x %x %x ", &type , &n0, &n1, 
            &n2, &c0, &c1 );
          FaceTypes->InsertValue(k,type);
          FaceNodes->InsertComponent(k,0,n0);
          FaceNodes->InsertComponent(k,1,n1);
          FaceNodes->InsertComponent(k,2,n2);
          FaceNodes->InsertComponent(k,3,0);
          FaceCells->InsertComponent(k,0,c0);
          FaceCells->InsertComponent(k,1,c1);
          }
        else
          {
          sscanf( buf, " %x %x %x %x %x %x %x ", &type, &n0 , &n1,
            &n2, &n3, &c0, &c1 );
          FaceTypes->InsertValue(k,type);
          FaceNodes->InsertComponent(k,0,n0);
          FaceNodes->InsertComponent(k,1,n1);
          FaceNodes->InsertComponent(k,2,n2);
          FaceNodes->InsertComponent(k,3,n3);
          FaceCells->InsertComponent(k,0,c0);
          FaceCells->InsertComponent(k,1,c1);
          }
        }
      else if ( et == 2)
        {
        sscanf( buf, " %x %x %x %x ", &n0 , &n1, &c0, &c1 );
        FaceTypes->InsertValue(k,2);
        FaceNodes->InsertComponent(k,0,n0);
        FaceNodes->InsertComponent(k,1,n1);
        FaceNodes->InsertComponent(k,2,0);
        FaceNodes->InsertComponent(k,3,0);
        FaceCells->InsertComponent(k,0,c0);
        FaceCells->InsertComponent(k,1,c1);
        }
      else if ( et == 3)
        {
        sscanf( buf, " %x %x %x %x %x ", &n0 , &n1, &n2, &c0, &c1 );
        FaceTypes->InsertValue(k,3);
        FaceNodes->InsertComponent(k,0,n0);
        FaceNodes->InsertComponent(k,1,n1);
        FaceNodes->InsertComponent(k,2,n2);
        FaceNodes->InsertComponent(k,3,0);
        FaceCells->InsertComponent(k,0,c0);
        FaceCells->InsertComponent(k,1,c1);
        }
      else
        {
        sscanf( buf, " %x %x %x %x %x %x ", &n0 , &n1, &n2, &n3, &c0, &c1 );
        FaceTypes->InsertValue(k,4);
        FaceNodes->InsertComponent(k,0,n0);
        FaceNodes->InsertComponent(k,1,n1);
        FaceNodes->InsertComponent(k,2,n2);
        FaceNodes->InsertComponent(k,3,n3);
        FaceCells->InsertComponent(k,0,c0);
        FaceCells->InsertComponent(k,1,c1);
        }
      j = GoToNextEOL(j) +1;
      }
    }
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  Points->InsertPoint(0, 0.0 , 0.0 , 0.0);

  if (zi == 0)
    {
    NumberOfNodes = li;
    }
  else
    {
    j = GoToNextLeftParen(j)+1;
    j = GoToNextEOL(j) +1;
    float x,y,z;
    for (int k = fi; k <= li; k++)
      {
      GetStringToNextRightParenOrEOL( j, buf );
      if ( nd == 2)
        {
        sscanf( buf, " %f %f ", &x , &y );
        Points->InsertPoint(k, x, y, 0.0);
        }
      else
        {
        sscanf( buf, " %f %f %f", &x , &y, &z );
        Points->InsertPoint(k, x, y, z);
        }

      j = GoToNextEOL(j) +1;
      }
    }
  j = GoToNextRightParen(j)+1;
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);

  j = GoToNextLeftParen(j)+1;
  j = GoToNextASCIIHexDigit(j);

  for (int k=face_id0;k<=face_id1;k++)
    {
    GetStringToNextRightParenOrEOL( j, buf );
    int pid0, pid1;
    sscanf( buf, " %x %x ", &pid0 , &pid1 );

    FaceParents->InsertComponent(k, 0, pid0);
    FaceParents->InsertComponent(k, 1, pid1);
    FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    NumberOfFaceParentChildren++;

    j = GoToNextEOL(j) +1;
    }

  if (face_id1 >= NumberOfFaceParents)
    {
    NumberOfFaceParents = face_id1;
    }

  return GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationASCII(int ix)
{
  // Face Information
  char buf[120];
  int j = ix + 1;	

  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );

  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);
  NCGFaceKidId->InsertValue(NumberOfNCGFaceHeaders, KidId);
  NCGFaceParentId->InsertValue(NumberOfNCGFaceHeaders, ParentId);
  NCGFaceNumberOfFaces->InsertValue(NumberOfNCGFaceHeaders, NumberOfFacesNCG);

  j = GoToNextLeftParen(j)+1;
  j = GoToNextASCIIHexDigit(j);

  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    GetStringToNextRightParenOrEOL( j, buf );
    int child, parent;
    sscanf( buf, " %x %x ", &child , &parent );
    NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    j = GoToNextEOL(j) +1;
    NumberOfNCGFaces++;
    }

  NumberOfNCGFaceHeaders++;
  return GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationASCII(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);
  NCGNodeZoneId->InsertValue(NumberOfNCGNodeHeaders, ZoneId);
  NCGNodeNumberOfNodesNCG->InsertValue(NumberOfNCGNodeHeaders,
    NumberOfNodesNCG);

  j = GoToNextLeftParen(j)+1;
  j = GoToNextASCIIHexDigit(j);

  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    GetStringToNextRightParenOrEOL( j, buf );
    float x,y,z;
    int NodeId;
    if (GridDimension == 3)
      {
      sscanf( buf, " %d %f %f %f ", &NodeId, &x , &y, &z );
      }
    else
      {
      sscanf( buf, " %d %f %f ", &NodeId, &x , &y );
      z = 0;
      }

    NCGNodeIds->InsertValue(NumberOfNCGNodes, NodeId);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 0, x);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 1, y);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 2, z);

    j = GoToNextEOL(j) +1;
    NumberOfNCGNodes++;
    }

  NumberOfNCGNodeHeaders++;
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodeFlagsASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::Command54(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetZoneSectionsASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = GoToNextLeftParen(j)+1;
  j = GoToNextASCIIHexDigit(j);

  int psf0, psf1;
  for (int k = fi; k <= li; k++)
    {
    GetStringToNextRightParenOrEOL( j, buf );
    sscanf( buf, " %x %x ", &psf0 , &psf1 );
    PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    j = GoToNextEOL(j) +1;
    }

  if ( li >= NumberOfPeriodicShadowFaces)
    {
    NumberOfPeriodicShadowFaces = li;
    }

  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetGridSizeASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPartitionASCII(int ix)
{
  return GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  CellTreeParentCellId0->InsertValue(NumberOfCellTrees, fid0);
  CellTreeParentCellId1->InsertValue(NumberOfCellTrees, fid1);
  CellTreeParentZoneId->InsertValue(NumberOfCellTrees, pzid);
  CellTreeChildZoneId->InsertValue(NumberOfCellTrees, czid);

  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetAsciiInteger(j);
    j = GoPastAsciiInteger(j);
    CellTreesNumberOfKids->InsertValue(NumberOfCellTreeParents, NumberOfKids);
    CellTreesKidsIndex->InsertValue(NumberOfCellTreeParents, 
      NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetAsciiInteger(j);
      j = GoPastAsciiInteger(j);
      CellTreesKids->InsertValue(NumberOfCellTreeKids, Kid);
      NumberOfCellTreeKids++;
      }
    NumberOfCellTreeParents++;
    }

  NumberOfCellTrees++;
  return GoToNextSectionASCII(j);
}

int vtkFLUENTReader::GetFaceTreeASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	

  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  FaceTreeParentFaceId0->InsertValue(NumberOfFaceTrees, fid0);
  FaceTreeParentFaceId1->InsertValue(NumberOfFaceTrees, fid1);
  FaceTreeParentZoneId->InsertValue(NumberOfFaceTrees, pzid);
  FaceTreeChildZoneId->InsertValue(NumberOfFaceTrees, czid);

  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetAsciiInteger(j);
    j = GoPastAsciiInteger(j);
    FaceTreesNumberOfKids->InsertValue(NumberOfFaceTreeParents, NumberOfKids);
    FaceTreesKidsIndex->InsertValue(NumberOfFaceTreeParents, 
      NumberOfFaceTreeKids);
    for(int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetAsciiInteger(j);
      j = GoPastAsciiInteger(j);
      FaceTreesKids->InsertValue(NumberOfFaceTreeKids, Kid);
      NumberOfFaceTreeKids++;
      }

    NumberOfFaceTreeParents++;
    }

  NumberOfFaceTrees++;
  return GoToNextSectionASCII(j);
}

//-----------------------------------------------------------------------------

int vtkFLUENTReader::GetVariablesSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2037)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCortexVariablesSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2038)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDomainVariablesSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2064)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if ( zi != 0)
    {
    CellZones->InsertValue(NumberOfCellZones, zi);
    NumberOfCellZones++;
    }

  if ( et != 0)
    {
    for (int i = fi; i <= li; i++)
      {
      CellTypes->InsertValue(i, et);
      }
    }
  else
    { // Mixed Cells
    j = GoToNextLeftParen(j)+1;
    for (int i = fi; i <= li; i++)
      {
      CellTypes->InsertValue(i, GetBinaryInteger(j));
      j = j + 4;
      }
    }

  j++;
  return GoToNextSectionSinglePrecision( j, "2012)");	
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  j = GoToNextLeftParen(j)+1;

  if ( et == 2)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2, 0);
      FaceNodes->InsertComponent(i,3, 0);
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if ( et == 3)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,3, 0);
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if ( et == 4)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,3,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else
    { // Mixed Faces
    for (int i = fi; i <= li; i++)
      {
      int ft = GetBinaryInteger(j);
      j = j + 4;
      FaceTypes->InsertValue(i, ft);
      if ( ft == 2)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2, 0);
        FaceNodes->InsertComponent(i,3, 0);
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 3)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,3, 0);
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 4)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,3,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      }
    }
  return GoToNextSectionSinglePrecision( j, "2013)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  Points->InsertPoint(0, 0.0 , 0.0 , 0.0);
  j = GoToNextLeftParen(j)+1;

  float x,y,z;
  for(int k = fi; k <= li; k++)
    {
    if ( nd == 2)
      {
      x = GetBinaryFloat(j);
      j = j+4;
      y = GetBinaryFloat(j);
      j = j+4;
      Points->InsertPoint(k, x, y, 0);
      }
    else
      {
      x = GetBinaryFloat(j);
      j = j+4;
      y = GetBinaryFloat(j);
      j = j+4;
      z = GetBinaryFloat(j);
      j = j+4;
      Points->InsertPoint(k, x, y, z);
      }
    }

  return GoToNextSectionSinglePrecision( j, "2010)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);
  j = GoToNextLeftParen(j)+1;

  int pid0, pid1;
  for (int k = face_id0; k <= face_id1; k++)
    {
    pid0 = GetBinaryInteger(j);
    j = j + 4;
    pid1 = GetBinaryInteger(j);
    j = j + 4;
    FaceParents->InsertComponent(k, 0, pid0);
    FaceParents->InsertComponent(k, 1, pid1);
    FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    NumberOfFaceParentChildren++;
    }

  if ( face_id1 >= NumberOfFaceParents)
    {
    NumberOfFaceParents = face_id1;
    }

  return GoToNextSectionSinglePrecision( j, "2061)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);
  NCGFaceKidId->InsertValue(NumberOfNCGFaceHeaders, KidId);
  NCGFaceParentId->InsertValue(NumberOfNCGFaceHeaders, ParentId);
  NCGFaceNumberOfFaces->InsertValue(NumberOfNCGFaceHeaders, 
    NumberOfFacesNCG);

  j = GoToNextLeftParen(j)+1;
  int child,parent;
  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    child = GetBinaryInteger(j);
    j = j + 4;
    parent = GetBinaryInteger(j);
    j = j + 4;
    NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    NumberOfNCGFaces++;
    }

  NumberOfNCGFaceHeaders++;
  return GoToNextSectionSinglePrecision( j, "2062)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationSinglePrecision(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);
  NCGNodeZoneId->InsertValue(NumberOfNCGNodeHeaders, ZoneId);
  NCGNodeNumberOfNodesNCG->InsertValue(NumberOfNCGNodeHeaders, 
    NumberOfNodesNCG);
  j = GoToNextLeftParen(j)+1;

  float x,y,z;
  int NodeId;
  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    NodeId = GetBinaryInteger(j);
    j = j + 4;
    x = GetBinaryFloat(j);
    j = j + 4;
    y = GetBinaryFloat(j);
    j = j + 4;
    if ( GridDimension == 3)
      {
      z = GetBinaryFloat(j);
      j = j + 4;
      }
    else
      {
      z = 0.0;
      }

    NCGNodeIds->InsertValue(NumberOfNCGNodes, NodeId);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 0, x);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 1, y);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 2, z);
    NumberOfNCGNodes++;
    }

  NumberOfNCGNodeHeaders++;
  return GoToNextSectionSinglePrecision( j, "2063)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodeFlagsSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2041)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetZoneSectionsSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2039)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = GoToNextLeftParen(j)+1;

  int psf0, psf1;
  for (int k = fi; k <= li; k++)
    {
    psf0 = GetBinaryInteger(j);
    j = j + 4;
    psf1 = GetBinaryInteger(j);
    j = j + 4;
    PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    }

  if ( li >= NumberOfPeriodicShadowFaces)
    {
    NumberOfPeriodicShadowFaces = li;
    }

  return GoToNextSectionSinglePrecision( j, "2018)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetGridSizeSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2033)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPartitionSinglePrecision(int ix)
{
  return GoToNextSectionSinglePrecision( ix, "2040)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  CellTreeParentCellId0->InsertValue(NumberOfCellTrees, fid0);
  CellTreeParentCellId1->InsertValue(NumberOfCellTrees, fid1);
  CellTreeParentZoneId->InsertValue(NumberOfCellTrees, pzid);
  CellTreeChildZoneId->InsertValue(NumberOfCellTrees, czid);
  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetBinaryInteger(j);
    j = j + 4;
    CellTreesNumberOfKids->InsertValue(NumberOfCellTreeParents, NumberOfKids);
    CellTreesKidsIndex->InsertValue(NumberOfCellTreeParents, 
      NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetBinaryInteger(j);
      j = j + 4;
      CellTreesKids->InsertValue(NumberOfCellTreeKids, Kid);
      NumberOfCellTreeKids++;
      }

    NumberOfCellTreeParents++;
    }

  NumberOfCellTrees++;
  return GoToNextSectionSinglePrecision( j, "2058)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceTreeSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  FaceTreeParentFaceId0->InsertValue(NumberOfFaceTrees, fid0);
  FaceTreeParentFaceId1->InsertValue(NumberOfFaceTrees, fid1);
  FaceTreeParentZoneId->InsertValue(NumberOfFaceTrees, pzid);
  FaceTreeChildZoneId->InsertValue(NumberOfFaceTrees, czid);
  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetBinaryInteger(j);
    j = j + 4;
    FaceTreesNumberOfKids->InsertValue(NumberOfFaceTreeParents, NumberOfKids);
    FaceTreesKidsIndex->InsertValue(NumberOfFaceTreeParents, 
      NumberOfFaceTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetBinaryInteger(j);
      j = j + 4;
      FaceTreesKids->InsertValue(NumberOfFaceTreeKids, Kid);
      NumberOfFaceTreeKids++;
      }
    NumberOfFaceTreeParents++;
    }
  NumberOfFaceTrees++;
  return GoToNextSectionSinglePrecision( j, "2059)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetVariablesDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3037)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCortexVariablesDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3038)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDomainVariablesDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3064)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  if ( zi != 0)
    {
    CellZones->InsertValue(NumberOfCellZones, zi);
    NumberOfCellZones++;
    }

  if ( et != 0)
    {
    for ( int i = fi; i <= li; i++)
      {
      CellTypes->InsertValue(i, et);
      }
    }
  else
    { // Mixed Cells
    j = GoToNextLeftParen(j)+1;
    for (int i = fi; i <= li; i++)
      {
      CellTypes->InsertValue(i, GetBinaryInteger(j));
      j = j + 4;
      }
    }

  j++;
  return GoToNextSectionDoublePrecision( j, "3012)");	
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  j = GoToNextLeftParen(j)+1;

  if (et == 2)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2, 0);
      FaceNodes->InsertComponent(i,3, 0);
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if (et == 3)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,3, 0);
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if (et == 4)
    {
    for (int i = fi; i <= li; i++)
      {
      FaceTypes->InsertValue(i, et);
      FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
      j = j + 4;
      FaceNodes->InsertComponent(i,3,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
      j = j + 4;
      FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else
    { // Mixed Faces
    for (int i = fi; i <= li; i++)
      {
      int ft = GetBinaryInteger(j);
      j = j + 4;
      FaceTypes->InsertValue(i, ft);
      if ( ft == 2)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2, 0);
        FaceNodes->InsertComponent(i,3, 0);
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 3)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,3, 0);
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 4)
        {
        FaceNodes->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,2,GetBinaryInteger(j));
        j = j + 4;
        FaceNodes->InsertComponent(i,3,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,0,GetBinaryInteger(j));
        j = j + 4;
        FaceCells->InsertComponent(i,1,GetBinaryInteger(j));
        j = j + 4;
        }
      }
    }
  return GoToNextSectionDoublePrecision( j, "3013)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  Points->InsertPoint(0, 0.0 , 0.0 , 0.0);
  j = GoToNextLeftParen(j)+1;

  float x,y,z;
  for (int k = fi; k <= li; k++)
    {
    if ( nd == 2)
      {
      x = GetBinaryDouble(j);
      j = j+8;
      y = GetBinaryDouble(j);
      j = j+8;
      Points->InsertPoint(k, x, y, 0);
      }
    else
      {
      x = GetBinaryDouble(j);
      j = j+8;
      y = GetBinaryDouble(j);
      j = j+8;
      z = GetBinaryDouble(j);
      j = j+8;
      Points->InsertPoint(k, x, y, z);
      }
    }
  return GoToNextSectionSinglePrecision( j, "3010)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);
  j = GoToNextLeftParen(j)+1;

  int pid0, pid1;
  for (int k = face_id0; k <= face_id1; k++)
    {
    pid0 = GetBinaryInteger(j);
    j = j + 4;
    pid1 = GetBinaryInteger(j);
    j = j + 4;
    FaceParents->InsertComponent(k, 0, pid0);
    FaceParents->InsertComponent(k, 1, pid1);
    FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    NumberOfFaceParentChildren++;
    }

  if (face_id1 >= NumberOfFaceParents)
    {
    NumberOfFaceParents = face_id1;
    }

  return GoToNextSectionDoublePrecision( j, "3061)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationDoublePrecision(int ix)
{
  // Face Information
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;

  GetStringToNextRightParen( j, buf );
  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);
  NCGFaceKidId->InsertValue(NumberOfNCGFaceHeaders, KidId);
  NCGFaceParentId->InsertValue(NumberOfNCGFaceHeaders, ParentId);
  NCGFaceNumberOfFaces->InsertValue(NumberOfNCGFaceHeaders, NumberOfFacesNCG);

  j = GoToNextLeftParen(j)+1;
  int child,parent;
  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    child = GetBinaryInteger(j);
    j = j + 4;
    parent = GetBinaryInteger(j);
    j = j + 4;
    NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    NumberOfNCGFaces++;
    }

  NumberOfNCGFaceHeaders++;
  return GoToNextSectionDoublePrecision( j, "3062)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationDoublePrecision(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);
  NCGNodeZoneId->InsertValue(NumberOfNCGNodeHeaders, ZoneId);
  NCGNodeNumberOfNodesNCG->InsertValue(NumberOfNCGNodeHeaders, 
    NumberOfNodesNCG);
  j = GoToNextLeftParen(j)+1;

  float x,y,z;
  int NodeId;
  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    NodeId = GetBinaryInteger(j);
    j = j + 4;
    x = GetBinaryDouble(j);
    j = j + 8;
    y = GetBinaryDouble(j);
    j = j + 8;
    if (GridDimension == 3)
      {
      z = GetBinaryDouble(j);
      j = j + 8;
      }
    else
      {
      z = 0.0;
      }

    NCGNodeIds->InsertValue(NumberOfNCGNodes, NodeId);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 0, x);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 1, y);
    NCGNodes->InsertComponent(NumberOfNCGNodes, 2, z);
    NumberOfNCGNodes++;
    }

  NumberOfNCGNodeHeaders++;
  return GoToNextSectionDoublePrecision( j, "3063)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodeFlagsDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3041)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetZoneSectionsDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3039)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );
  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = GoToNextLeftParen(j)+1;

  int psf0, psf1;
  for (int k = fi;k <= li; k++)
    {
    psf0 = GetBinaryInteger(j);
    j = j + 4;
    psf1 = GetBinaryInteger(j);
    j = j + 4;
    PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    }

  if ( li >= NumberOfPeriodicShadowFaces)
    {
    NumberOfPeriodicShadowFaces = li;
    }

  return GoToNextSectionDoublePrecision( j, "3018)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetGridSizeDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3033)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPartitionDoublePrecision(int ix)
{
  return GoToNextSectionDoublePrecision( ix, "3040)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);
  CellTreeParentCellId0->InsertValue(NumberOfCellTrees, fid0);
  CellTreeParentCellId1->InsertValue(NumberOfCellTrees, fid1);
  CellTreeParentZoneId->InsertValue(NumberOfCellTrees, pzid);
  CellTreeChildZoneId->InsertValue(NumberOfCellTrees, czid);

  j = GoToNextLeftParen(j)+1;
  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetBinaryInteger(j);
    j = j + 4;
    CellTreesNumberOfKids->InsertValue(NumberOfCellTreeParents, NumberOfKids);
    CellTreesKidsIndex->InsertValue(NumberOfCellTreeParents, 
      NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetBinaryInteger(j);
      j = j + 4;
      CellTreesKids->InsertValue(NumberOfCellTreeKids, Kid);
      NumberOfCellTreeKids++;
      }
    NumberOfCellTreeParents++;
    }
  NumberOfCellTrees++;
  return GoToNextSectionDoublePrecision( j, "3058)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceTreeDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParen(j)+1;
  GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);
  FaceTreeParentFaceId0->InsertValue(NumberOfFaceTrees, fid0);
  FaceTreeParentFaceId1->InsertValue(NumberOfFaceTrees, fid1);
  FaceTreeParentZoneId->InsertValue(NumberOfFaceTrees, pzid);
  FaceTreeChildZoneId->InsertValue(NumberOfFaceTrees, czid);
  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = GetBinaryInteger(j);
    j = j + 4;
    FaceTreesNumberOfKids->InsertValue(NumberOfFaceTreeParents, NumberOfKids);
    FaceTreesKidsIndex->InsertValue(NumberOfFaceTreeParents, 
      NumberOfFaceTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = GetBinaryInteger(j);
      j = j + 4;
      FaceTreesKids->InsertValue(NumberOfFaceTreeKids, Kid);
      NumberOfFaceTreeKids++;
      }
    NumberOfFaceTreeParents++;
    }
  NumberOfFaceTrees++;
  return GoToNextSectionDoublePrecision( j, "3059)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetGridDimension(int ix)
{
  char b2[2];
  b2[0] = CaseFileBuffer[ix+3];
  b2[1] = 0;
  GridDimension = atoi(b2);

  return GoToNextRightParen(ix);
}


//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataComment(int ix)
{
  return GoToNextRightParenData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataHeader(int ix)
{
  return GoToNextRightParenData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataGridDimension(int ix)
{
  char b2[2];
  b2[0] = DataFileBuffer[ix+3];
  b2[1] = 0;
  GridDimension = atoi(b2);

  return GoToNextRightParenData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataMachineConfiguration(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataGridSizeASCII(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataVariablesASCII(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownASCII313(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextRightParenData(int ix)
{
  int i = ix;
  while ( DataFileBuffer[i] != ')' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextLeftParenData(int ix)
{
  int i = ix;
  while ( DataFileBuffer[i] != '(' )
    {
    i++;
    }
    return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionASCIIData(int ix)
{
  int i = ix + 1;
  int level = 0;

  while ( !((level == 0) && (DataFileBuffer[i] == ')')))
    {
    if ( DataFileBuffer[i] == '(')
      {
      level++;
      }
    if (DataFileBuffer[i] == ')')
      {
      level--;
      }
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionSinglePrecisionData(int ix, char buf[])
{
  int i = ix + 1;
  while ( !((DataFileBuffer[i] == buf[0]) 
    && (DataFileBuffer[i+1] == buf[1]) 
    && (DataFileBuffer[i+2] == buf[2]) 
    && (DataFileBuffer[i+3] == buf[3]) 
    && (DataFileBuffer[i+4] == buf[4])))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionDoublePrecisionData(int ix, char buf[])
{
  int i = ix + 1;
  while( !((DataFileBuffer[i] == buf[0]) 
    && (DataFileBuffer[i+1] == buf[1]) 
    && (DataFileBuffer[i+2] == buf[2]) 
    && (DataFileBuffer[i+3] == buf[3]) 
    && (DataFileBuffer[i+4] == buf[4])))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataASCII(int ix)
{
  char buf[120];
  int j = ix + 1;	
  float x;
  j = GoToNextLeftParenData(j)+1;
  GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = GoToNextLeftParenData(j)+1;

  if (DataPass == 1)
    {
    if (IsCellZoneId(zid))
      {
      if (IsNewVariable(ssid))
        {
        VariableIds->InsertValue(NumberOfVariables, ssid);
        VariableSizes->InsertValue(NumberOfVariables, size);
        NumberOfVariables++;
        }
      }
    }
  else 
    {
    if (IsCellZoneId(zid))
      {
      int index = GetVariableIndex(ssid);
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          GetStringToNextRightParenOrEOLData( j, buf );
          sscanf( buf, " %f ", &x );
          j = GoToNextEOLData(j) +1;
          CellData[index]->InsertComponent( i-1, k, x); 
          }
        }
      }
    }
  return GoToNextSectionASCIIData(j);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParenData(j)+1;
  GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = GoToNextLeftParenData(j)+1;

  if ( DataPass == 1)
    {
    if ( IsCellZoneId(zid))
      {
      if ( IsNewVariable(ssid))
        {
        VariableIds->InsertValue(NumberOfVariables, ssid);
        VariableSizes->InsertValue(NumberOfVariables, size);
        NumberOfVariables++;
        }
      }
    }
  else
    {
    if ( IsCellZoneId(zid))
      {
      int index = GetVariableIndex(ssid);
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          CellData[index]->InsertComponent( i-1, k, GetBinaryFloatData(j)); 
          j = j + 4;
          }
        }
      }
    }
  return GoToNextSectionSinglePrecisionData( j, "2300)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;	
  j = GoToNextLeftParenData(j)+1;
  GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = GoToNextLeftParenData(j)+1;

  if ( DataPass == 1)
    {
    if ( IsCellZoneId(zid))
      {
      if ( IsNewVariable(ssid))
        {
        VariableIds->InsertValue(NumberOfVariables, ssid);
        VariableSizes->InsertValue(NumberOfVariables, size);
        NumberOfVariables++;
        }
      }
    }
  else
    {
    if ( IsCellZoneId(zid))
      {
      int index = GetVariableIndex(ssid);
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          CellData[index]->InsertComponent( i-1, k, GetBinaryDoubleData(j)); 
          j = j + 8;
          }
        }
      }
    }
  return GoToNextSectionSinglePrecisionData( j, "3300)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownSinglePrecision2301(int ix)
{
  return GoToNextSectionSinglePrecisionData( ix, "2301)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownSinglePrecision2302(int ix)
{
  return GoToNextSectionSinglePrecisionData( ix, "2302)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownSinglePrecision2313(int ix)
{
  return GoToNextSectionSinglePrecisionData( ix, "2313)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownDoublePrecision3301(int ix)
{
  return GoToNextSectionDoublePrecisionData( ix, "3301)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownDoublePrecision3302(int ix)
{
  return GoToNextSectionDoublePrecisionData( ix, "3302)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownDoublePrecision3313(int ix)
{
  return GoToNextSectionDoublePrecisionData( ix, "3313)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownASCII301(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownASCII302(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetUnknownASCII303(int ix)
{
  return GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataUnknownASCII(int ix)
{
  int j = ix + 1;	
  j = GoToNextLeftParenData(j)+1;
  j = GoToNextLeftParenData(j)+1;
  j = GoToNextRightParenData(j)+1;
  j = GoToNextRightParenData(j)+1;
  return j;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetStringToNextRightParenData(int ix, char buf[] )
{
  // Copy contents between ( ) into buffer
  int j = ix;
  int k=0;
  while ( !(DataFileBuffer[j] == ')'))
    {
    buf[k] = DataFileBuffer[j];
    j++;
    k++;
    }
  buf[k] = 0;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::IsCellZoneId(int zi)
{
  int flag = 0;
  for (int i = 0; i < NumberOfCellZones; i++)
    {
    if ( zi == CellZones->GetValue(i))
      {
      flag = 1;
      }
    }
  return flag;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::IsNewVariable(int ssid)
{
  int flag = 1;
  for (int i = 0; i < NumberOfVariables; i++)
    {
    if ( ssid == VariableIds->GetValue(i))
      {
      flag = 0;
      }
    }
  return flag;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetVariableIndex(int ssid)
{
  int index = 0;
  for (int i = 0; i < NumberOfVariables; i++)
    {
    if ( ssid == VariableIds->GetValue(i))
      {
      index = i;
      }
    }
  return index;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetBinaryIntegerData(int ix)
{
  union mix_i
    {
    int i;
    char c[4];
    }mi= {1};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = DataFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
float vtkFLUENTReader::GetBinaryFloatData(int ix)
{
  union mix_i
    {
    float i;
    char c[4];
    } mi = {1.0};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = DataFileBuffer[ix+k];
      }
    }	
  return mi.i;
}

//-----------------------------------------------------------------------------
double vtkFLUENTReader::GetBinaryDoubleData(int ix)
{
  union mix_i
    {
    double i;
    char c[8];
    } mi= {1.0};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 7; k >= 0; k--)
      {
      mi.c[k] = DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 7; k++)
      {
      mi.c[7-k] = DataFileBuffer[ix+k];
      }
    }	
  return mi.i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::IsASCIICharacterHexDigit(int ix)
{
  if ( (CaseFileBuffer[ix] >= 0x30) && (CaseFileBuffer[ix] <= 0x39))
    {
    return 1;
    }
  else if ( (CaseFileBuffer[ix] >= 0x41) && (CaseFileBuffer[ix] <= 0x46))
    {
    return 1;
    }
  else if ( (CaseFileBuffer[ix] >= 0x61) && (CaseFileBuffer[ix] <= 0x66))
    {
    return 1;
    }
  else
    {
    return 0;
    }
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextASCIIHexDigit(int ix)
{
  int i = ix;
  while (! IsASCIICharacterHexDigit(i))
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextRightParen(int ix)
{
  int i = ix;
  while (CaseFileBuffer[i] != ')' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextLeftParen(int ix)
{
  int i = ix;
  while (CaseFileBuffer[i] != '(' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextEOL(int ix)
{
  int i = ix;
  while (CaseFileBuffer[i] != 0x0a )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionASCII(int ix)
{
  int i = ix + 1;
  int level = 0;
  while ( !((level == 0) && (CaseFileBuffer[i] == ')')))
    {
    if (CaseFileBuffer[i] == '(')
      {
      level++;
      }
    if (CaseFileBuffer[i] == ')')
      {
      level--;
      }
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetStringToNextRightParen(int ix, char buf[] )
{
  // Copy contents between ( ) into buffer
  int j = ix;
  int k=0;
  while ( !(CaseFileBuffer[j] == ')'))
    {
    buf[k] = CaseFileBuffer[j];
    j++;
    k++;
    }
  buf[k] = 0;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetStringToNextRightParenOrEOL(int ix, char buf[] )
{
  // Copy contents between ( ) into buffer
  int j = ix;
  int k=0;
  while ( !((CaseFileBuffer[j] == 0x0a)||(CaseFileBuffer[j] == ')')))
    {
    buf[k] = CaseFileBuffer[j];
    j++;
    k++;
    }
  buf[k] = 0;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetMixedCellTypes(int ix, int fi, int li)
{
  int j = ix;
  char c[2];
  for (int i = fi; i <= li; i++)
    {
    j = GoToNextASCIIHexDigit(j);
    cout << "i = " << i << ", et = " << CaseFileBuffer[j] << endl;
    c[0] = CaseFileBuffer[j];
    c[1] = 0;
    CellTypes->InsertValue(i, atoi(c));
    j++;
    }
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionSinglePrecision(int ix, char buf[])
{
  int i = ix + 1;
  while( !((CaseFileBuffer[i] == buf[0]) 
    && (CaseFileBuffer[i+1] == buf[1]) 
    && (CaseFileBuffer[i+2] == buf[2]) 
    && (CaseFileBuffer[i+3] == buf[3]) 
    && (CaseFileBuffer[i+4] == buf[4]) ))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionDoublePrecision(int ix, char buf[])
{
  int i = ix + 1;
  while ( !((CaseFileBuffer[i] == buf[0]) 
    && (CaseFileBuffer[i+1] == buf[1]) 
    && (CaseFileBuffer[i+2] == buf[2]) 
    && (CaseFileBuffer[i+3] == buf[3]) 
    && (CaseFileBuffer[i+4] == buf[4])))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetBinaryInteger(int ix)
{
  union mix_i
    {
    int i;
    char c[4];
    } mi = {1};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = CaseFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
float vtkFLUENTReader::GetBinaryFloat(int ix)
{
  union mix_i
    {
    float i;
    char c[4];
    } mi = {1.0};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = CaseFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
double vtkFLUENTReader::GetBinaryDouble(int ix)
{
  union mix_i
    {
    double i;
    char c[8];
    } mi = {1.0};

  if ( LittleEndianFlag == 1)
    {
    for (int k = 7; k >= 0; k--)
      {
      mi.c[k] = CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 7; k++)
      {
      mi.c[7-k] = CaseFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetAsciiInteger(int ix)
{
  int j = ix;
  int first = GoToNextASCIIHexDigit(j);
  j = first;
  while ( IsASCIICharacterHexDigit(++j));
  int second = j-1;
  int nchar = second-first+1;
  char *buf;
  buf = new char[nchar+1];
  buf[nchar] = 0;
  j = first;

  for (int i = 0; i < nchar; i++)
    {
    buf[i] = CaseFileBuffer[j];
    j++;
    }
  return atoi(buf);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoPastAsciiInteger(int ix)
{
  int j = ix;
  int first = GoToNextASCIIHexDigit(j);
  j = first;
  while(IsASCIICharacterHexDigit(++j));
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextEOLData(int ix)
{
  int i = ix;
  while(DataFileBuffer[i] != 0x0a )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetStringToNextRightParenOrEOLData(int ix, char buf[] )
{
  // Copy contents between ( ) into buffer
  int j = ix;
  int k=0;
  while ( !((DataFileBuffer[j] == 0x0a)||(DataFileBuffer[j] == ')'))) 
    {
    buf[k] = DataFileBuffer[j];
    j++;
    k++;
    }
  buf[k] = 0;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::CreateVTKObjects(void)
{
  this->NumberOfCellFields = 0;
  this->NumberOfCellComponents = 0;
  this->FileStream = NULL;
  this->NumberOfCells = 0;
  this->CellDataInfo = NULL;
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->SetNumberOfInputPorts(0);

  this->CaseFileBuffer = NULL;
  this->DataFileBuffer = NULL;
  this->CaseFileBufferLength = 0;
  this->DataFileBufferLength = 0;
  this->GridDimension = 0;
  this->NumberOfNodes = 0;
  this->NumberOfFaces = 0;
  this->NumberOfFaceParents = 0;
  this->NumberOfPeriodicShadowFaces = 0;
  this->NumberOfCellZones = 0;
  this->NumberOfVariables = 0;
  this->LittleEndianFlag = 1;
  this->NumberOfFaceTrees = 0;
  this->NumberOfFaceTreeKids = 0;
  this->NumberOfFaceTreeParents = 0;
  this->LastFaceTreeParent = 0;
  this->NumberOfCellTrees = 0;
  this->NumberOfCellTreeKids = 0;
  this->NumberOfCellTreeParents = 0;
  this->LastCellTreeParent = 0;
  this->NumberOfNCGFaceHeaders = 0;
  this->NumberOfNCGFaces = 0;
  this->NumberOfNCGNodeHeaders = 0;
  this->NumberOfNCGNodes = 0;
  this->DataPass = 0;
  this->NumberOfFaceParentChildren = 0;

  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->Points = vtkPoints::New();
  this->CellTypes = vtkIntArray::New();
  this->CellFaces = vtkIntArray::New();
  this->CellFacesClean = vtkIntArray::New();
  this->CellFacesClean->SetNumberOfComponents(6);
  this->FaceTypes = vtkIntArray::New();
  this->FaceNodes = vtkIntArray::New();
  this->FaceNodes->SetNumberOfComponents(4);
  this->FaceCells = vtkIntArray::New();
  this->FaceCells->SetNumberOfComponents(2);
  this->FaceParents = vtkIntArray::New();
  this->FaceParents->SetNumberOfComponents(2);
  this->PeriodicShadowFaces = vtkIntArray::New();
  this->PeriodicShadowFaces->SetNumberOfComponents(2);
  this->FaceTreesNumberOfKids = vtkIntArray::New();
  this->FaceTreesKids = vtkIntArray::New();
  this->FaceTreesKidsIndex = vtkIntArray::New();
  this->CellTreesNumberOfKids = vtkIntArray::New();
  this->CellTreesKids = vtkIntArray::New();
  this->CellTreesKidsIndex = vtkIntArray::New();
  this->FaceTreeParentFaceId0 = vtkIntArray::New();
  this->FaceTreeParentFaceId1 = vtkIntArray::New();
  this->FaceTreeParentZoneId = vtkIntArray::New();
  this->FaceTreeChildZoneId = vtkIntArray::New();
  this->FaceTreeParentTable = vtkIntArray::New();
  this->CellTreeParentCellId0 = vtkIntArray::New();
  this->CellTreeParentCellId1 = vtkIntArray::New();
  this->CellTreeParentZoneId = vtkIntArray::New();
  this->CellTreeChildZoneId = vtkIntArray::New();
  this->CellTreeParentTable = vtkIntArray::New();
  this->NCGFaceKidId = vtkIntArray::New();
  this->NCGFaceParentId = vtkIntArray::New();
  this->NCGFaceNumberOfFaces = vtkIntArray::New();
  this->NCGFaceChild = vtkIntArray::New();
  this->NCGFaceParent = vtkIntArray::New();
  this->NCGNodeZoneId = vtkIntArray::New();
  this->NCGNodeNumberOfNodesNCG = vtkIntArray::New();
  this->NCGNodes = vtkDoubleArray::New();
  this->NCGNodes->SetNumberOfComponents(2);
  this->NCGNodeIds = vtkIntArray::New();
  this->CellNumberOfFaces = vtkIntArray::New();
  this->FaceKidFlags = vtkIntArray::New();
  this->FaceParentFlags = vtkIntArray::New();
  this->CellIndex = vtkIntArray::New();
  this->InterfaceFaceChildFlags = vtkIntArray::New();
  this->FaceParentsChildren = vtkIntArray::New();
  this->NCGFaceChildFlags = vtkIntArray::New();
  this->CellParentFlags = vtkIntArray::New();
  this->ATriangle = vtkTriangle::New();
  this->AQuad = vtkQuad::New();
  this->ATetra = vtkTetra::New();
  this->APyramid = vtkPyramid::New();
  this->AWedge = vtkWedge::New();
  this->AHexahedron = vtkHexahedron::New();
  this->CellZones = vtkIntArray::New();
  this->VariableIds = vtkIntArray::New();
  this->VariableSizes = vtkIntArray::New();
  this->CellData = NULL;
  this->Mesh = vtkUnstructuredGrid::New();

  ObjectsFlag = 1;
}


//-----------------------------------------------------------------------------
void vtkFLUENTReader::DeleteVTKObjects(void)
{
  delete [] CaseFileBuffer;
  delete [] DataFileBuffer;
  Points->Delete();
  CellTypes->Delete();
  CellFaces->Delete();
  CellFacesClean->Delete();

  FaceTypes->Delete();
  FaceNodes->Delete();
  FaceCells->Delete();
  FaceParents->Delete();
  PeriodicShadowFaces->Delete();
  FaceTreesNumberOfKids->Delete();
  FaceTreesKids->Delete();
  FaceTreesKidsIndex->Delete();
  CellTreesNumberOfKids->Delete();
  CellTreesKids->Delete();
  CellTreesKidsIndex->Delete();
  FaceTreeParentFaceId0->Delete();
  FaceTreeParentFaceId1->Delete();
  FaceTreeParentZoneId->Delete();
  FaceTreeChildZoneId->Delete();
  FaceTreeParentTable->Delete();
  CellTreeParentCellId0->Delete();
  CellTreeParentCellId1->Delete();
  CellTreeParentZoneId->Delete();
  CellTreeChildZoneId->Delete();
  CellTreeParentTable->Delete();

  NCGFaceKidId->Delete();
  NCGFaceParentId->Delete();
  NCGFaceNumberOfFaces->Delete();
  NCGFaceChild->Delete();
  NCGFaceParent->Delete();
  NCGNodeZoneId->Delete();
  NCGNodeNumberOfNodesNCG->Delete();
  NCGNodes->Delete();
  NCGNodeIds->Delete();
  CellNumberOfFaces->Delete();
  FaceKidFlags->Delete();
  FaceParentFlags->Delete();
  CellIndex->Delete();
  InterfaceFaceChildFlags->Delete();
  FaceParentsChildren->Delete();
  NCGFaceChildFlags->Delete();
  CellParentFlags->Delete();
  ATriangle->Delete();
  AQuad->Delete();
  ATetra->Delete();
  APyramid->Delete();
  AWedge->Delete();
  AHexahedron->Delete();
  CellZones->Delete();
  VariableIds->Delete();
  VariableSizes->Delete();
  ObjectsFlag = 0;
}
