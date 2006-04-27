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

// mccdo
  for ( int i = 0; i < this->NumberOfVariables; i++ ) 
    {
    if ( this->CellData[ i ] )
      {
      this->CellData[ i ]->Delete();
      }
    }
// mccdo
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
  this->Mesh->SetPoints(this->Points);   //mccdo
  this->Points->Delete();   //mccdo
  this->MakeFaceTreeParentTable();
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

  for ( int i=0; i < this->NumberOfVariables; i++ ) 
    {
    int variableId = this->VariableIds->GetValue(i);
    int numberOfComponents = this->VariableSizes->GetValue(i);
    this->CellData[ i ] = vtkDoubleArray::New();
    this->CellData[ i ]->SetName(VariableNames[variableId]);
    this->CellData[ i ]->SetNumberOfComponents(numberOfComponents);
    }

  this->DataPass = 2;
  this->ParseDataFile();  // Getting Data

/*mccdo
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
*/
  this->NumberOfCellArrays = this->NumberOfCellFields;
  //mccdo this->Mesh->SetPoints(this->Points);
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
/* mccdo
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
*/
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
  //mccdo
  this->FileStream->seekg(0, ios::end); 
  this->CaseFileBufferLength = this->FileStream->tellg();
  this->FileStream->seekg(0, ios::beg);
  this->CaseFileBuffer = new char[this->CaseFileBufferLength];
  this->FileStream->read(this->CaseFileBuffer, this->CaseFileBufferLength);
  this->FileStream->close();
  //mccdo

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
  delete [] CaseFileBuffer; // mccdo
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::MakeFaceTreeParentTable(void)
{
/* mccdo
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
*/
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadFaceParentFlags(void)
{
/* mccdo  // Initialize
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

  this->FaceTreeParentFaceId0->Delete();
  this->FaceTreeParentFaceId1->Delete();
*/
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadInterfaceFaceChildFlags(void)
{
/* mccdo
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
*/
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::LoadNCGFaceChildFlags(void)
{
/* mccdo
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
  this->NCGFaceChild->Delete();
*/
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

  for(int i=1;i<=this->NumberOfCells;i++)
    {
    for (int j = 0; j < 6; j++)
      {
      try
        {
        face[j] = CellFaces[ i ].at(j); // mccdo (int)this->CellFacesClean->GetComponent(i, j);
        }
      catch ( ... )
        {
        face[j] = 0;
        }
      }

    for (int j = 0; j < 6; j++)
      {
      if ( (face[j]!=0) && ((int)this->FaceCells->GetComponent(face[j],0) == i))
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

    if(this->CellTypes->GetValue(i) == 1)
      {
      int cnt = 0;
      for (int j = 0; j < 3; j++)
        {
        for (int k = 0; k < 2; k++)
          {
          tempNode[cnt] = (int)this->FaceNodes->GetComponent(face[j], k);
          cnt++;
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

      for (int j = 0; j < 3; j++)
        {
        this->ATriangle->GetPointIds()->SetId( j, node[j]);
        }

      if(this->CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->ATriangle->GetCellType(), 
          this->ATriangle->GetPointIds());
        }
      }
    else if(this->CellTypes->GetValue(i) == 2)
      {
      //*************************************
      //   Tetrahedral Cell Type
      //*************************************
      int cnt = 0;
      for(int j = 0; j < 4; j++)
        {
        for(int k = 0; k < 3; k++)
          {
          tempNode[cnt] = (int)this->FaceNodes->GetComponent(face[j], k);
          cnt++;
          }
        }

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

      if ( (tempNode[3]!=node[0]) && (tempNode[3]!=node[1]) 
        && (tempNode[3]!=node[2])) 
        {
        node[3] = tempNode[3];
        }
      else if ( (tempNode[4]!=node[0]) && (tempNode[4]!=node[1])
        && (tempNode[4]!=node[2])) 
        {
        node[3] = tempNode[4];
        }
      else if ( (tempNode[5]!=node[0]) && (tempNode[5]!=node[1])
        && (tempNode[5]!=node[2])) 
        {
        node[3] = tempNode[5];
        }
      else if ( (tempNode[6]!=node[0]) && (tempNode[6]!=node[1])
        && (tempNode[6]!=node[2]))
        {
        node[3] = tempNode[6];
        }
      else if ( (tempNode[7]!=node[0]) && (tempNode[7]!=node[1])
        && (tempNode[7]!=node[2])) 
        {
        node[3] = tempNode[7];
        }
      else if ( (tempNode[8]!=node[0]) && (tempNode[8]!=node[1])
        && (tempNode[8]!=node[2])) 
        {
        node[3] = tempNode[8];
        }
      else if ( (tempNode[9]!=node[0]) && (tempNode[9]!=node[1])
        && (tempNode[9]!=node[2]))
        {
        node[3] = tempNode[9];
        }
      else if ( (tempNode[10]!=node[0]) && (tempNode[10]!=node[1])
        && (tempNode[10]!=node[2])) 
        {
        node[3] = tempNode[10];
        }
      else
        {
        node[3] = tempNode[11];
        }

      for (int j = 0; j < 4; j++)
        {
        this->ATetra->GetPointIds()->SetId( j, node[j]);
        }

      if (CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->ATetra->GetCellType(), 
          this->ATetra->GetPointIds());
        }
      } 
    else if (this->CellTypes->GetValue(i) == 3)
      {
      //*************************************
      //   Quadrilateral Cell Type
      //*************************************

      int cnt = 0;
      for(int j = 0; j < 4; j++)
        {
        for(int k = 0; k < 2; k++)
          {
          tempNode[cnt] = (int)this->FaceNodes->GetComponent(face[j], k);
          cnt++;
          }
        }
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

      if ( (tempNode[2]!=node[0]) && (tempNode[2]!=node[1])
        && (tempNode[3]!=node[0]) && (tempNode[3]!=node[1]))
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

      if ( (tempNode[4]!=node[0]) && (tempNode[4]!=node[1])
        && (tempNode[5]!=node[0]) && (tempNode[5]!=node[1]))
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

      if ( (tempNode[6]!=node[0]) && (tempNode[6]!=node[1])
        && (tempNode[7]!=node[0]) && (tempNode[7]!=node[1]))
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

      for (int j = 0; j < 4; j++)
        {
        this->AQuad->GetPointIds()->SetId( j, node[j]);
        }

      if (CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->AQuad->GetCellType(), 
          this->AQuad->GetPointIds());
        }
      }
    else if (this->CellTypes->GetValue(i) == 4)
      {
      //*************************************
      //   Hexahedral Cell Type
      //*************************************
      int cnt = 0;
      for(int j = 0; j < 6; j++)
        {
        for(int k = 0; k < 4; k++)
          {
          tempNode[cnt] = (int)this->FaceNodes->GetComponent(face[j], k);
          cnt++;
          }
        }

      if (spinFace[0] > 0)
        {
        for(int j = 0; j < 4; j++)
          {
          node[j] = tempNode[j];
          }
        }
      else
        {
        for(int j = 0; j < 4; j++)
          {
          node[j] = tempNode[3-j];
          }
        }

      int frontFaceNode[4];
      int topFaceNode[4];
      int rightFaceNode[4];
      int leftFaceNode[4];
      int bottomFaceNode[4];

      for (int j = 0; j < 4; j++)
        {
        frontFaceNode[j] = 0;
        topFaceNode[j] = 0;
        rightFaceNode[j] = 0;
        leftFaceNode[j] = 0;
        bottomFaceNode[j] = 0;
        }


      if (((tempNode[4]==node[0])||(tempNode[5]==node[0])
        ||(tempNode[6]==node[0])||(tempNode[7]==node[0]))
        && ((tempNode[4]==node[1])||(tempNode[5]==node[1])
        ||(tempNode[6]==node[1])||(tempNode[7]==node[1])))
        {
        for (int j = 0; j < 4; j++)
          {
          rightFaceNode[j] = tempNode[j+4];
          }
        }
      else if (((tempNode[4]==node[0])||(tempNode[5]==node[0])
        ||(tempNode[6]==node[0])||(tempNode[7]==node[0]))
        && ((tempNode[4]==node[3])||(tempNode[5]==node[3])
        ||(tempNode[6]==node[3])||(tempNode[7]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          frontFaceNode[j] = tempNode[j+4];
          }
        }
      else if (((tempNode[4]==node[2])||(tempNode[5]==node[2])
        ||(tempNode[6]==node[2])||(tempNode[7]==node[2])) 
        && ((tempNode[4]==node[3])||(tempNode[5]==node[3])
        ||(tempNode[6]==node[3])||(tempNode[7]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          leftFaceNode[j] = tempNode[j+4];
          }
        }
      else if (((tempNode[4]==node[1])||(tempNode[5]==node[1])
        ||(tempNode[6]==node[1])||(tempNode[7]==node[1])) 
        && ((tempNode[4]==node[2])||(tempNode[5]==node[2])
        ||(tempNode[6]==node[2])||(tempNode[7]==node[2])))
        {
        for (int j = 0; j < 4; j++)
          {
          bottomFaceNode[j] = tempNode[j+4];
          }
        }
      else
        {
        for (int j = 0; j < 4; j++)
          {
          topFaceNode[j] = tempNode[j+4];
          }
        }

      if (((tempNode[8]==node[0])||(tempNode[9]==node[0])
        ||(tempNode[10]==node[0])||(tempNode[11]==node[0]))
        && ((tempNode[8]==node[1])||(tempNode[9]==node[1])
        ||(tempNode[10]==node[1])||(tempNode[11]==node[1])))
        {
        for (int j = 0; j < 4; j++)
          {
          rightFaceNode[j] = tempNode[j+8];
          }
        }
      else if (((tempNode[8]==node[0])||(tempNode[9]==node[0])
        ||(tempNode[10]==node[0])||(tempNode[11]==node[0]))
        && ((tempNode[8]==node[3])||(tempNode[9]==node[3])
        ||(tempNode[10]==node[3])||(tempNode[11]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          frontFaceNode[j] = tempNode[j+8];
          }
        }
      else if (((tempNode[8]==node[2])||(tempNode[9]==node[2])
        ||(tempNode[10]==node[2])||(tempNode[11]==node[2]))
        && ((tempNode[8]==node[3])||(tempNode[9]==node[3])
        ||(tempNode[10]==node[3])||(tempNode[11]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          leftFaceNode[j] = tempNode[j+8];
          }
        }
      else if (((tempNode[8]==node[1])||(tempNode[9]==node[1])
        ||(tempNode[10]==node[1])||(tempNode[11]==node[1]))
        && ((tempNode[8]==node[2])||(tempNode[9]==node[2])
        ||(tempNode[10]==node[2])||(tempNode[11]==node[2])))
        {
        for (int j = 0; j < 4; j++)
          {
          bottomFaceNode[j] = tempNode[j+8];
          }
        }
      else
        {
        for (int j = 0; j < 4; j++)
          {
          topFaceNode[j] = tempNode[j+8];
          }
        }

      if (((tempNode[12]==node[0])||(tempNode[13]==node[0])
        ||(tempNode[14]==node[0])||(tempNode[15]==node[0]))
        && ((tempNode[12]==node[1])||(tempNode[13]==node[1])
        ||(tempNode[14]==node[1])||(tempNode[15]==node[1])))
        {
        for (int j = 0; j < 4; j++)
          {
          rightFaceNode[j] = tempNode[j+12];
          }
        }
      else if (((tempNode[12]==node[0])||(tempNode[13]==node[0])
        ||(tempNode[14]==node[0])||(tempNode[15]==node[0]))
        && ((tempNode[12]==node[3])||(tempNode[13]==node[3])
        ||(tempNode[14]==node[3])||(tempNode[15]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          frontFaceNode[j] = tempNode[j+12];
          }
        }
      else if (((tempNode[12]==node[2])||(tempNode[13]==node[2])
        ||(tempNode[14]==node[2])||(tempNode[15]==node[2]))
        && ((tempNode[12]==node[3])||(tempNode[13]==node[3])
        ||(tempNode[14]==node[3])||(tempNode[15]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          leftFaceNode[j] = tempNode[j+12];
          }
        }
      else if (((tempNode[12]==node[1])||(tempNode[13]==node[1])
        ||(tempNode[14]==node[1])||(tempNode[15]==node[1]))
        && ((tempNode[12]==node[2])||(tempNode[13]==node[2])
        ||(tempNode[14]==node[2])||(tempNode[15]==node[2])))
        {
        for (int j = 0; j < 4; j++)
          {
          bottomFaceNode[j] = tempNode[j+12];
          }
        }
      else
        {
        for (int j = 0; j < 4; j++)
          {
          topFaceNode[j] = tempNode[j+12];
          }
        }

      if (((tempNode[16]==node[0])||(tempNode[17]==node[0])
        ||(tempNode[18]==node[0])||(tempNode[19]==node[0]))
        && ((tempNode[16]==node[1])||(tempNode[17]==node[1])
        ||(tempNode[18]==node[1])||(tempNode[19]==node[1])))
        {
        for (int j = 0; j < 4; j++)
          {
          rightFaceNode[j] = tempNode[j+16];
          }
        }
      else if (((tempNode[16]==node[0])||(tempNode[17]==node[0])
        ||(tempNode[18]==node[0])||(tempNode[19]==node[0]))
        && ((tempNode[16]==node[3])||(tempNode[17]==node[3])
        ||(tempNode[18]==node[3])||(tempNode[19]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          frontFaceNode[j] = tempNode[j+16];
          }
        }
      else if (((tempNode[16]==node[2])||(tempNode[17]==node[2])
        ||(tempNode[18]==node[2])||(tempNode[19]==node[2]))
        && ((tempNode[16]==node[3])||(tempNode[17]==node[3])
        ||(tempNode[18]==node[3])||(tempNode[19]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          leftFaceNode[j] = tempNode[j+16];
          }
        }
      else if (((tempNode[16]==node[1])||(tempNode[17]==node[1])
        ||(tempNode[18]==node[1])||(tempNode[19]==node[1]))
        && ((tempNode[16]==node[2])||(tempNode[17]==node[2])
        ||(tempNode[18]==node[2])||(tempNode[19]==node[2])))
        {
        for (int j = 0; j < 4; j++)
          {
          bottomFaceNode[j] = tempNode[j+16];
          }
        }
      else
        {
        for (int j = 0; j < 4; j++)
          {
          topFaceNode[j] = tempNode[j+16];
          }
        }

      if (((tempNode[20]==node[0])||(tempNode[21]==node[0])
        ||(tempNode[22]==node[0])||(tempNode[23]==node[0]))
        && ((tempNode[20]==node[1])||(tempNode[21]==node[1])
        ||(tempNode[22]==node[1])||(tempNode[23]==node[1])))
        {
        for (int j = 0; j < 4; j++)
          {
          rightFaceNode[j] = tempNode[j+20];
          }
        }
      else if (((tempNode[20]==node[0])||(tempNode[21]==node[0])
        ||(tempNode[22]==node[0])||(tempNode[23]==node[0]))
        && ((tempNode[20]==node[3])||(tempNode[21]==node[3])
        ||(tempNode[22]==node[3])||(tempNode[23]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          frontFaceNode[j] = tempNode[j+20];
          }
        }
      else if (((tempNode[20]==node[2])||(tempNode[21]==node[2])
        ||(tempNode[22]==node[2])||(tempNode[23]==node[2])) 
        && ((tempNode[20]==node[3])||(tempNode[21]==node[3])
        ||(tempNode[22]==node[3])||(tempNode[23]==node[3])))
        {
        for (int j = 0; j < 4; j++)
          {
          leftFaceNode[j] = tempNode[j+20];
          }
        }
      else if (((tempNode[20]==node[1])||(tempNode[21]==node[1])
        ||(tempNode[22]==node[1])||(tempNode[23]==node[1])) 
        && ((tempNode[20]==node[2])||(tempNode[21]==node[2])
        ||(tempNode[22]==node[2])||(tempNode[23]==node[2])))
        {
        for (int j = 0; j < 4; j++)
          {
          bottomFaceNode[j] = tempNode[j+20];
          }
        }
      else
        {
        for (int j = 0; j < 4; j++)
          {
          topFaceNode[j] = tempNode[j+20];
          }
        }

      if (((topFaceNode[0]==rightFaceNode[0])
        ||(topFaceNode[0]==rightFaceNode[1])
        ||(topFaceNode[0]==rightFaceNode[2])
        ||(topFaceNode[0]==rightFaceNode[3]))
        &&((topFaceNode[0]==frontFaceNode[0])
        ||(topFaceNode[0]==frontFaceNode[1])
        ||(topFaceNode[0]==frontFaceNode[2])
        ||(topFaceNode[0]==frontFaceNode[3])))
        {
        node[4] = topFaceNode[0];
        }
      else if (((topFaceNode[0]==rightFaceNode[0])
        ||(topFaceNode[0]==rightFaceNode[1])
        ||(topFaceNode[0]==rightFaceNode[2])
        ||(topFaceNode[0]==rightFaceNode[3]))
        &&((topFaceNode[0]==bottomFaceNode[0])
        ||(topFaceNode[0]==bottomFaceNode[1])
        ||(topFaceNode[0]==bottomFaceNode[2])
        ||(topFaceNode[0]==bottomFaceNode[3])))
        {
        node[5] = topFaceNode[0];
        }
      else if (((topFaceNode[0]==leftFaceNode[0])
        ||(topFaceNode[0]==leftFaceNode[1])
        ||(topFaceNode[0]==leftFaceNode[2])
        ||(topFaceNode[0]==leftFaceNode[3]))
        &&((topFaceNode[0]==bottomFaceNode[0])
        ||(topFaceNode[0]==bottomFaceNode[1])
        ||(topFaceNode[0]==bottomFaceNode[2])
        ||(topFaceNode[0]==bottomFaceNode[3])))
        {
        node[6] = topFaceNode[0];
        }
      else 
        {
        node[7] = topFaceNode[0];
        }

      if (((topFaceNode[1]==rightFaceNode[0])
        ||(topFaceNode[1]==rightFaceNode[1])
        ||(topFaceNode[1]==rightFaceNode[2])
        ||(topFaceNode[1]==rightFaceNode[3]))
        &&((topFaceNode[1]==frontFaceNode[0])
        ||(topFaceNode[1]==frontFaceNode[1])
        ||(topFaceNode[1]==frontFaceNode[2])
        ||(topFaceNode[1]==frontFaceNode[3])))
        {
        node[4] = topFaceNode[1];
        }
      else if (((topFaceNode[1]==rightFaceNode[0])
        ||(topFaceNode[1]==rightFaceNode[1])
        ||(topFaceNode[1]==rightFaceNode[2])
        ||(topFaceNode[1]==rightFaceNode[3]))
        &&((topFaceNode[1]==bottomFaceNode[0])
        ||(topFaceNode[1]==bottomFaceNode[1])
        ||(topFaceNode[1]==bottomFaceNode[2])
        ||(topFaceNode[1]==bottomFaceNode[3])))
        {
        node[5] = topFaceNode[1];
        }
      else if (((topFaceNode[1]==leftFaceNode[0])
        ||(topFaceNode[1]==leftFaceNode[1])
        ||(topFaceNode[1]==leftFaceNode[2])
        ||(topFaceNode[1]==leftFaceNode[3]))
        &&((topFaceNode[1]==bottomFaceNode[0])
        ||(topFaceNode[1]==bottomFaceNode[1])
        ||(topFaceNode[1]==bottomFaceNode[2])
        ||(topFaceNode[1]==bottomFaceNode[3])))
        {
        node[6] = topFaceNode[1];
        }
      else
        {
        node[7] = topFaceNode[1];
        }

      if (((topFaceNode[2]==rightFaceNode[0])
        ||(topFaceNode[2]==rightFaceNode[1])
        ||(topFaceNode[2]==rightFaceNode[2])
        ||(topFaceNode[2]==rightFaceNode[3]))
        &&((topFaceNode[2]==frontFaceNode[0])
        ||(topFaceNode[2]==frontFaceNode[1])
        ||(topFaceNode[2]==frontFaceNode[2])
        ||(topFaceNode[2]==frontFaceNode[3])))
        {
        node[4] = topFaceNode[2];
        }
      else if (((topFaceNode[2]==rightFaceNode[0])
        ||(topFaceNode[2]==rightFaceNode[1])
        ||(topFaceNode[2]==rightFaceNode[2])
        ||(topFaceNode[2]==rightFaceNode[3]))
        &&((topFaceNode[2]==bottomFaceNode[0])
        ||(topFaceNode[2]==bottomFaceNode[1])
        ||(topFaceNode[2]==bottomFaceNode[2])
        ||(topFaceNode[2]==bottomFaceNode[3])))
        {
        node[5] = topFaceNode[2];
        }
      else if (((topFaceNode[2]==leftFaceNode[0])
        ||(topFaceNode[2]==leftFaceNode[1])
        ||(topFaceNode[2]==leftFaceNode[2])
        ||(topFaceNode[2]==leftFaceNode[3]))
        &&((topFaceNode[2]==bottomFaceNode[0])
        ||(topFaceNode[2]==bottomFaceNode[1])
        ||(topFaceNode[2]==bottomFaceNode[2])
        ||(topFaceNode[2]==bottomFaceNode[3])))
        {
        node[6] = topFaceNode[2];
        }
      else
        {
        node[7] = topFaceNode[2];
        }

      if (((topFaceNode[3]==rightFaceNode[0])
        ||(topFaceNode[3]==rightFaceNode[1])
        ||(topFaceNode[3]==rightFaceNode[2])
        ||(topFaceNode[3]==rightFaceNode[3]))
        &&((topFaceNode[3]==frontFaceNode[0])
        ||(topFaceNode[3]==frontFaceNode[1])
        ||(topFaceNode[3]==frontFaceNode[2])
        ||(topFaceNode[3]==frontFaceNode[3])))
        {
        node[4] = topFaceNode[3];
        }
      else if (((topFaceNode[3]==rightFaceNode[0])
        ||(topFaceNode[3]==rightFaceNode[1])
        ||(topFaceNode[3]==rightFaceNode[2])
        ||(topFaceNode[3]==rightFaceNode[3])) 
        &&((topFaceNode[3]==bottomFaceNode[0])
        ||(topFaceNode[3]==bottomFaceNode[1])
        ||(topFaceNode[3]==bottomFaceNode[2])
        ||(topFaceNode[3]==bottomFaceNode[3])))
        {
        node[5] = topFaceNode[3];
        }
      else if (((topFaceNode[3]==leftFaceNode[0])
        ||(topFaceNode[3]==leftFaceNode[1])
        ||(topFaceNode[3]==leftFaceNode[2])
        ||(topFaceNode[3]==leftFaceNode[3]))
        &&((topFaceNode[3]==bottomFaceNode[0])
        ||(topFaceNode[3]==bottomFaceNode[1])
        ||(topFaceNode[3]==bottomFaceNode[2])
        ||(topFaceNode[3]==bottomFaceNode[3])))
        {
        node[6] = topFaceNode[3];
        }
      else
        {
        node[7] = topFaceNode[3];
        }

      for (int j = 0; j < 8; j++)
        {
        this->AHexahedron->GetPointIds()->SetId( j, node[j]);
        }

      if (CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->AHexahedron->GetCellType(), 
          this->AHexahedron->GetPointIds());
        }
      }
    else if (CellTypes->GetValue(i) == 5)
      {
      //*************************************
      //   Pyramid Cell Type
      //*************************************
      if (this->FaceTypes->GetValue(face[0]) == 4)
        {
        if (spinFace[0] > 0)
          {
          for (int j = 0; j < 4; j++)
            {
            node[j] = (int)this->FaceNodes->GetComponent(face[0], j);
            }
          }
        else
          {
          for (int j = 0; j < 4; j++)
            {
            node[3-j] = (int)this->FaceNodes->GetComponent(face[0], j);
            }
          }
        for (int j = 0; j < 3; j++)
          {
          tempNode[j] = (int)this->FaceNodes->GetComponent(face[1], j);
          }
        }
      else if (this->FaceTypes->GetValue(face[1]) == 4)
        {
        if (spinFace[1] > 0)
          {
          for (int j = 0; j < 4; j++)
            {
            node[j] = (int)this->FaceNodes->GetComponent(face[1], j);
            }
          }
        else
          {
          for (int j = 0; j < 4; j++)
            {
            node[3-j] = (int)this->FaceNodes->GetComponent(face[1], j);
            }
          }
        for (int j = 0; j < 3; j++)
          {
          tempNode[j] = (int)this->FaceNodes->GetComponent(face[0], j);
          }
        }
      else if (this->FaceTypes->GetValue(face[2]) == 4)
        {
        if (spinFace[2] > 0)
          {
          for (int j = 0; j < 4; j++)
            {
            node[j] = (int)this->FaceNodes->GetComponent(face[2], j);
            }
          }
        else
          {
          for (int j = 0; j < 4; j++)
            {
            node[3-j] = (int)this->FaceNodes->GetComponent(face[2], j);
            }
          }
        for (int j = 0; j < 3; j++)
          {
          tempNode[j] = (int)this->FaceNodes->GetComponent(face[0], j);
          }
        }
      else if (this->FaceTypes->GetValue(face[3]) == 4)
        {
        if (spinFace[3] > 0)
          {
          for (int j = 0; j < 4; j++)
            {
            node[j] = (int)this->FaceNodes->GetComponent(face[3], j);
            }
          }
        else
          {
          for (int j = 0; j < 4; j++)
            {
            node[3-j] = (int)this->FaceNodes->GetComponent(face[3], j);
            }
          }
        for (int j = 0; j < 3; j++)
          {
          tempNode[j] = (int)this->FaceNodes->GetComponent(face[0], j);
          }
        }
      else
        {
        if (spinFace[4] > 0)
          {
          for (int j = 0; j < 4; j++)
            {
            node[j] = (int)this->FaceNodes->GetComponent(face[4], j);
            }
          }
        else
          {
          for (int j = 0; j < 4; j++)
            {
            node[3-j] = (int)this->FaceNodes->GetComponent(face[4], j);
            }
          }
        for (int j = 0; j < 3; j++)
          {
          tempNode[j] = (int)this->FaceNodes->GetComponent(face[0], j);
          }
        }

      if ((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])
        &&(tempNode[0]!=node[2])&&(tempNode[0]!=node[3]))
        {
        node[4] = tempNode[0];
        }
      else if ((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])
        &&(tempNode[1]!=node[2])&&(tempNode[1]!=node[3]))
        {
        node[4] = tempNode[1];
        }
      else
        {
        node[4] = tempNode[2];
        }

      for (int j = 0; j < 5; j++)
        {
        this->APyramid->GetPointIds()->SetId( j, node[j]);
        }

      if (this->CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->APyramid->GetCellType(),
          this->APyramid->GetPointIds());
        }

      }
    else if (this->CellTypes->GetValue(i) == 6)
      {
      //*************************************
      //   Wedge Cell Type
      //*************************************
      if (this->FaceTypes->GetValue(face[0]) == 4)
        {
        if (spinFace[0] > 0)
          {
          node[0] = (int)this->FaceNodes->GetComponent(face[0], 0);
          node[1] = (int)this->FaceNodes->GetComponent(face[0], 1);
          node[4] = (int)this->FaceNodes->GetComponent(face[0], 2);
          node[3] = (int)this->FaceNodes->GetComponent(face[0], 3);
          }
        else
          {
          node[3] = (int)this->FaceNodes->GetComponent(face[0], 0);
          node[4] = (int)this->FaceNodes->GetComponent(face[0], 1);
          node[1] = (int)this->FaceNodes->GetComponent(face[0], 2);
          node[0] = (int)this->FaceNodes->GetComponent(face[0], 3);
          }
        }
      else if (this->FaceTypes->GetValue(face[1]) == 4)
        {
        if (spinFace[1] > 0)
          {
          node[0] = (int)this->FaceNodes->GetComponent(face[1], 0);
          node[1] = (int)this->FaceNodes->GetComponent(face[1], 1);
          node[4] = (int)this->FaceNodes->GetComponent(face[1], 2);
          node[3] = (int)this->FaceNodes->GetComponent(face[1], 3);
          }
        else
          {
          node[3] = (int)this->FaceNodes->GetComponent(face[1], 0);
          node[4] = (int)this->FaceNodes->GetComponent(face[1], 1);
          node[1] = (int)this->FaceNodes->GetComponent(face[1], 2);
          node[0] = (int)this->FaceNodes->GetComponent(face[1], 3);
          }
        }
      else if  (this->FaceTypes->GetValue(face[2]) == 4)
        {
        if (spinFace[2] > 0) 
          {
          node[0] = (int)this->FaceNodes->GetComponent(face[2], 0);
          node[1] = (int)this->FaceNodes->GetComponent(face[2], 1);
          node[4] = (int)this->FaceNodes->GetComponent(face[2], 2);
          node[3] = (int)this->FaceNodes->GetComponent(face[2], 3);
          }
        else
          {
          node[3] = (int)this->FaceNodes->GetComponent(face[2], 0);
          node[4] = (int)this->FaceNodes->GetComponent(face[2], 1);
          node[1] = (int)this->FaceNodes->GetComponent(face[2], 2);
          node[0] = (int)this->FaceNodes->GetComponent(face[2], 3);
          }
        }
      else if (this->FaceTypes->GetValue(face[3]) == 4)
        {
        if (spinFace[3] > 0)
          {
          node[0] = (int)this->FaceNodes->GetComponent(face[3], 0);
          node[1] = (int)this->FaceNodes->GetComponent(face[3], 1);
          node[4] = (int)this->FaceNodes->GetComponent(face[3], 2);
          node[3] = (int)this->FaceNodes->GetComponent(face[3], 3);
          }
        else
          {
          node[3] = (int)this->FaceNodes->GetComponent(face[3], 0);
          node[4] = (int)this->FaceNodes->GetComponent(face[3], 1);
          node[1] = (int)this->FaceNodes->GetComponent(face[3], 2);
          node[0] = (int)this->FaceNodes->GetComponent(face[3], 3);
          }
        }
      else
        {
        if (spinFace[4] > 0)
          {
          node[0] = (int)this->FaceNodes->GetComponent(face[4], 0);
          node[1] = (int)this->FaceNodes->GetComponent(face[4], 1);
          node[4] = (int)this->FaceNodes->GetComponent(face[4], 2);
          node[3] = (int)this->FaceNodes->GetComponent(face[4], 3);
          }
        else
          {
          node[3] = (int)this->FaceNodes->GetComponent(face[4], 0);
          node[4] = (int)this->FaceNodes->GetComponent(face[4], 1);
          node[1] = (int)this->FaceNodes->GetComponent(face[4], 2);
          node[0] = (int)this->FaceNodes->GetComponent(face[4], 3);
          }
        }

      int trf[2];
      int index = 0;
      if ( this->FaceTypes->GetValue(face[0]) == 3)
        {
        trf[index] = face[0];
        index++;
        }
      if ( this->FaceTypes->GetValue(face[1]) == 3)
        {
        trf[index] = face[1];
        index++;
        }
      if ( this->FaceTypes->GetValue(face[2]) == 3)
        {
        trf[index] = face[2];
        index++;
        }
      if ( this->FaceTypes->GetValue(face[3]) == 3)
        {
        trf[index] = face[3];
        index++;
        }
      if ( this->FaceTypes->GetValue(face[4]) == 3)
        {
        trf[index] = face[4];
        index++;
        }

      int cnt = 0;
      for (int j = 0; j < 2; j++)
        {
        for (int k = 0; k < 3; k++)
          {
          tempNode[cnt] = (int)this->FaceNodes->GetComponent(trf[j], k);
          cnt++;
          }
        }

      if (((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])
        &&(tempNode[0]!=node[4])&&(tempNode[0]!=node[3]))
        && ((tempNode[1]==node[0])||(tempNode[1]==node[1])) )
        {
        node[2] = tempNode[0];
        }
      else if (((tempNode[0]!=node[0])&&(tempNode[0]!=node[1])
        &&(tempNode[0]!=node[4])&&(tempNode[0]!=node[3]))
        && ((tempNode[1]==node[4])||(tempNode[1]==node[3])))
        {
        node[5] = tempNode[0];
        }

      if (((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])
        &&(tempNode[1]!=node[4])&&(tempNode[1]!=node[3]))
        && ((tempNode[0]==node[0])||(tempNode[0]==node[1])))
        {
        node[2] = tempNode[1];
        }
      else if (((tempNode[1]!=node[0])&&(tempNode[1]!=node[1])
        &&(tempNode[1]!=node[4])&&(tempNode[1]!=node[3]))
        && ((tempNode[0]==node[4])||(tempNode[0]==node[3])))
        {
        node[5] = tempNode[1];
        }

      if (((tempNode[2]!=node[0])&&(tempNode[2]!=node[1])
        &&(tempNode[2]!=node[4])&&(tempNode[2]!=node[3])) 
        && ((tempNode[1]==node[0])||(tempNode[1]==node[1]))) 
        {
        node[2] = tempNode[2];
        }
      else if (((tempNode[2]!=node[0])&&(tempNode[2]!=node[1])
        &&(tempNode[2]!=node[4])&&(tempNode[2]!=node[3]))
        && ((tempNode[1]==node[4])||(tempNode[1]==node[3])))
        {
        node[5] = tempNode[2];
        }

      if (((tempNode[3]!=node[0])&&(tempNode[3]!=node[1])
        &&(tempNode[3]!=node[4])&&(tempNode[3]!=node[3]))
        && ((tempNode[4]==node[0])||(tempNode[4]==node[1])))
        {
        node[2] = tempNode[3];
        }
      else if (((tempNode[3]!=node[0])&&(tempNode[3]!=node[1])
        &&(tempNode[3]!=node[4])&&(tempNode[3]!=node[3]))
        && ((tempNode[4]==node[4])||(tempNode[4]==node[3])))
        {
        node[5] = tempNode[3];
        }

      if (((tempNode[4]!=node[0])&&(tempNode[4]!=node[1])
        &&(tempNode[4]!=node[4])&&(tempNode[4]!=node[3])) 
        && ((tempNode[3]==node[0])||(tempNode[3]==node[1])))
        {
        node[2] = tempNode[4];
        }
      else if (((tempNode[4]!=node[0])&&(tempNode[4]!=node[1])
        &&(tempNode[4]!=node[4])&&(tempNode[4]!=node[3])) 
        && ((tempNode[3]==node[4])||(tempNode[3]==node[3])))
        {
        node[5] = tempNode[4];
        }

      if (((tempNode[5]!=node[0])&&(tempNode[5]!=node[1])
        &&(tempNode[5]!=node[4])&&(tempNode[5]!=node[3]))
        && ((tempNode[4]==node[0])||(tempNode[4]==node[1])))
        {
        node[2] = tempNode[5];
        }
      else if (((tempNode[5]!=node[0])&&(tempNode[5]!=node[1])
        &&(tempNode[5]!=node[4])&&(tempNode[5]!=node[3]))
        && ((tempNode[4]==node[4])||(tempNode[4]==node[3])))
        {
        node[5] = tempNode[5];
        }
      for (int j = 0; j < 6; j++)
        {
        this->AWedge->GetPointIds()->SetId( j, node[j]);
        }

      if (this->CellParentFlags.at(i) != true ) // mccdo
        {
        this->Mesh->InsertNextCell(this->AWedge->GetCellType(), 
          this->AWedge->GetPointIds());
        }
      }
    }
  // mccdo
  CellTypes->Delete();
  CellFaces.clear();
  CellParentFlags.clear();
  FaceTypes->Delete();
  FaceNodes->Delete();
  FaceCells->Delete();
  // mccdo
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellParentFlags(void)
{
/* mccdo
  // Initialize Array
  for (int i = 1; i <= this->NumberOfCells; i++)
    {
    this->CellParentFlags->InsertValue(i,0);
    }

  for (int i = 0; i < this->NumberOfCellTrees; i++)
    {
    for (int j = this->CellTreeParentCellId0->GetValue(i); 
      j <= this->CellTreeParentCellId1->GetValue(i); j++)
      {
      this->CellParentFlags->InsertValue(j,1);
      }
    }
*/
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellNumberOfFaces(void)
{
/* mccdo
  for (int i = 0; i <= this->NumberOfCells; i++)
    {
    this->CellNumberOfFaces->InsertValue( i, 0);
    }

  for (int i = 1; i <= this->NumberOfFaces; i++)
    {
    int c0 = (int)this->FaceCells->GetComponent( i, 0);
    int c1 = (int)this->FaceCells->GetComponent( i, 1);
    int nc0 = this->CellNumberOfFaces->GetValue(c0);
    int nc1 = this->CellNumberOfFaces->GetValue(c1);

    if ( c0 != 0)
      {
      nc0++;
      this->CellNumberOfFaces->InsertValue( c0, nc0);
      }

    if ( c1 != 0)
      {
      nc1++;
      this->CellNumberOfFaces->InsertValue( c1, nc1);
      }
    }
*/
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::LoadCellFaces(void)
{
  // Make an index array to determine where each cell is in the cell 
  // face array and ...
  // Make a temporary number of faces/cell array to keep track of 
  // where to put the faces within each block.
/*
  int index = 0;
  int *NumberOfFacesInCell;
  NumberOfFacesInCell = new int[this->NumberOfCells+1];

  for (int i = 1; i <= this->NumberOfCells; i++)
    {
    this->CellIndex->InsertValue(i, index);
    index = index + this->CellNumberOfFaces->GetValue(i);
    NumberOfFacesInCell[i] = 0;
    }

  this->CellIndex->InsertValue(0, 0);
  NumberOfFacesInCell[0] = 0;

  for (int i = 1; i <= this->NumberOfFaces; i++)
    {
    int c0 = (int)this->FaceCells->GetComponent(i,0);
    int c1 = (int)this->FaceCells->GetComponent(i,1);
    int nc0 = NumberOfFacesInCell[c0];
    int nc1 = NumberOfFacesInCell[c1];
    int ic0 = this->CellIndex->GetValue(c0);
    int ic1 = this->CellIndex->GetValue(c1);

    if ( c0 != 0)
      {
      this->CellFaces->InsertValue( ic0+nc0, i);
      nc0++;
      NumberOfFacesInCell[c0] = nc0;
      }

    if ( c1 != 0 )
      {
      this->CellFaces->InsertValue( ic1 + nc1, i);
      nc1++;
      NumberOfFacesInCell[c1] = nc1;
      }
    }

  delete [] NumberOfFacesInCell;
*/
   // mccdo
  for (int i = 1; i <= this->NumberOfFaces; i++)
    {
    int c0 = (int)this->FaceCells->GetComponent(i,0);
    int c1 = (int)this->FaceCells->GetComponent(i,1);
    if ( c0 != 0)
      {
      vtkstd::map< int, std::vector< int > >::iterator cellMap;
      cellMap = CellFaces.find( c0 );
      if ( cellMap == CellFaces.end() )
        {
        std::vector< int > temp;
        temp.push_back( i );
        CellFaces[ c0 ] = temp;
        }
      else
        {
        cellMap->second.push_back( i );
        }
      }

    if ( c1 != 0 )
      {
      vtkstd::map< int, std::vector< int > >::iterator cellMap;
      cellMap = CellFaces.find( c1 );
      if ( cellMap == CellFaces.end() )
        {
        std::vector< int > temp;
        temp.push_back( i );
        CellFaces[ c1 ] = temp;
        }
      else
        {
        cellMap->second.push_back( i );
        }
      }
    }
   // mccdo
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::RemoveExtraFaces(void)
{
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
  /* mccdo
  for (int i = 0; i <= this->NumberOfCells; i++)
    {
    for(int j = 0; j < 6; j++)
      {
      this->CellFacesClean->InsertComponent( i, j, 0);
      }
    }
  */

  for (int i = 1; i <= this->NumberOfCells; i++)
    {
    numberOfBadKids = 0;
    int cellType = this->CellTypes->GetValue(i);
    int numberOfFaces = CellFaces[ i ].size(); // mccdo this->CellNumberOfFaces->GetValue(i);

    if ( numberOfFaces > actualFaces[cellType])
      {
      // mccdo int ic = this->CellIndex->GetValue(i);
      for (int j = 0; j < numberOfFaces; j++)
        {
        // mccdo int face = this->CellFaces->GetValue(ic+j);
        // mccdo int parentFlag = this->FaceParentFlags->GetValue(face);
        // mccdo int ifChildFlag = this->InterfaceFaceChildFlags->GetValue(face);
        // mccdo
        int face = this->CellFaces[ i ].at( j );
        int parentFlag = this->FaceParentFlags[ face ];
        int ifChildFlag = InterfaceFaceChildFlags.at( face );
        // mccdo

        // mccdo int ncgFaceChildFlag = this->NCGFaceChildFlags->GetValue(face);
        // mccdo
        vtkstd::set< int >::iterator ncgChildIter;
        ncgChildIter = this->NCGFaceChildFlags.find( face );
        int ncgFaceChildFlag = 0;
        if ( ncgChildIter == NCGFaceChildFlags.end() )
          {
          ncgFaceChildFlag = 0;
          }
        else
          {
          ncgFaceChildFlag = 1;
          }
        // mccdo

        if (parentFlag == 1)
          {
          int startKid = 
            this->FaceTreesKidsIndex->GetValue(
            this->FaceTreeParentTable->GetValue(face));
          int endKid = this->FaceTreesKidsIndex->GetValue(
            this->FaceTreeParentTable->GetValue(face))
            + this->FaceTreesNumberOfKids->GetValue(
            this->FaceTreeParentTable->GetValue(face));
          for (int k = startKid; k < endKid; k++)
            {
            badKids[numberOfBadKids] = this->FaceTreesKids->GetValue(k);
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

      if ( (numberOfBadKids + actualFaces[cellType]) != numberOfFaces)
        {
        cout << " Problem in Face Reduction !!!! " << endl;
        cout << " Cell = " << i << endl;
        cout << " Problem - Number of Faces = " 
          << numberOfFaces << ", Actual Faces " 
          << actualFaces[CellTypes->GetValue(i)] 
          << ", Cell Type = " << CellTypes->GetValue(i) << endl;
        }

      // mccdo int idx = 0;
      std::vector< int >::iterator faceIter;
      faceIter =  CellFaces[ i ].begin();
      for (int j = 0; j < (int)CellFaces[ i ].size(); )
        {
        int bk = 0;
        // mccdo int face = CellFaces->GetValue( ic + j);
        int face = CellFaces[ i ].at( j );
        for (int m = 0; m < numberOfBadKids; m++)
          {
          if ( badKids[m] == face)
            {
            bk = 1;
            //remove this face from the vector
            faceIter = CellFaces[ i ].erase( faceIter );
            break;
            }
          }

        if (bk == 0)
          {
          // mccdo faces[idx] = face;
          // mccdo idx++;
          faceIter++;
          j++;
          }
        }
      }
    else
      {
      /* mccdo
      int idx = 0;
      int ic = this->CellIndex->GetValue(i);
      for (int j = 0; j < numberOfFaces; j++)
        {
        int face = this->CellFaces->GetValue(ic+j);
        faces[idx] = face;
        idx++;
        }
      */
      }
    /*
    for (int j = 0; j < actualFaces[cellType]; j++)
      {
      this->CellFacesClean->InsertComponent( i, j, faces[j]);
      }
    */
    }
  // mccdo
  this->FaceTreesNumberOfKids->Delete();
  this->FaceTreesKidsIndex->Delete();
  this->FaceTreeParentTable->Delete();
  this->FaceTreesKids->Delete();
  this->InterfaceFaceChildFlags.clear();
  this->FaceParentFlags.clear();
  this->NCGFaceChildFlags.clear();
  // mccdo
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::ParseDataFile(void)
{
  // mccdo
  this->DataFileStream->seekg(0, ios::end);
  this->DataFileBufferLength = this->DataFileStream->tellg();
  this->DataFileStream->seekg(0, ios::beg);
  this->DataFileBuffer = new char[this->DataFileBufferLength];
  this->DataFileStream->read(this->DataFileBuffer, this->DataFileBufferLength);
  // mccdo this->DataFileStream->close();
  // mccdo

  int bufptr = 0;
  while ( bufptr < this->DataFileBufferLength)
    {
    if ( this->DataFileBuffer[bufptr] == '(')
      {
      int ix = this->GetDataIndex(bufptr);
      bufptr = this->ExecuteDataTask(ix, bufptr);
      }
    bufptr++;
    }
  delete [] DataFileBuffer; // mccdo
  return;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::InitializeVariableNames ( void )
{
  this->VariableNames[1] = "PRESSURE";
  this->VariableNames[2] = "MOMENTUM";
  this->VariableNames[3] = "TEMPERATURE";
  this->VariableNames[4] = "ENTHALPY";
  this->VariableNames[5] = "TKE";
  this->VariableNames[6] = "TED";
  this->VariableNames[7] = "SPECIES";
  this->VariableNames[8] = "G";
  this->VariableNames[9] = "WSWIRL";
  this->VariableNames[10] = "DPMS_MASS";
  this->VariableNames[11] = "DPMS_MOM";
  this->VariableNames[12] = "DPMS_ENERGY";
  this->VariableNames[13] = "DPMS_SPECIES";
  this->VariableNames[14] = "DVOLUME_DT";
  this->VariableNames[15] = "BODY_FORCES";
  this->VariableNames[16] = "FMEAN";
  this->VariableNames[17] = "FVAR";
  this->VariableNames[18] = "MASS_FLUX";
  this->VariableNames[19] = "WALL_SHEAR";
  this->VariableNames[20] = "BOUNDARY_HEAT_FLUX";
  this->VariableNames[21] = "BOUNDARY_RAD_HEAT_FLUX";
  this->VariableNames[22] = "OLD_PRESSURE";
  this->VariableNames[23] = "POLLUT";
  this->VariableNames[24] = "DPMS_P1_S";
  this->VariableNames[25] = "DPMS_P1_AP";
  this->VariableNames[26] = "WALL_GAS_TEMPERATURE";
  this->VariableNames[27] = "DPMS_P1_DIFF";
  this->VariableNames[28] = "DR_SURF";
  this->VariableNames[29] = "W_M1";
  this->VariableNames[30] = "W_M2";
  this->VariableNames[31] = "DPMS_BURNOUT";
  this->VariableNames[32] = "DPMS_CONCENTRATION";
  this->VariableNames[33] = "PDF_MW";
  this->VariableNames[34] = "DPMS_WSWIRL";
  this->VariableNames[35] = "YPLUS";
  this->VariableNames[36] = "YPLUS_UTAU";
  this->VariableNames[37] = "WALL_SHEAR_SWIRL";
  this->VariableNames[38] = "WALL_T_INNER";
  this->VariableNames[39] = "POLLUT0";
  this->VariableNames[40] = "POLLUT1";
  this->VariableNames[41] = "WALL_G_INNER";
  this->VariableNames[42] = "PREMIXC";
  this->VariableNames[43] = "PREMIXC_T";
  this->VariableNames[44] = "PREMIXC_RATE";
  this->VariableNames[45] = "POLLUT2";
  this->VariableNames[46] = "POLLUT3";
  this->VariableNames[47] = "MASS_FLUX_M1";
  this->VariableNames[48] = "MASS_FLUX_M2";
  this->VariableNames[49] = "GRID_FLUX";
  this->VariableNames[50] = "DO_I";
  this->VariableNames[51] = "DO_RECON_I";
  this->VariableNames[52] = "DO_ENERGY_SOURCE";
  this->VariableNames[53] = "DO_IRRAD";
  this->VariableNames[54] = "DO_QMINUS";
  this->VariableNames[55] = "DO_IRRAD_OLD";
  this->VariableNames[56] = "DO_IWX";
  this->VariableNames[57] = "DO_IWY";
  this->VariableNames[58] = "DO_IWZ";
  this->VariableNames[59] = "MACH";
  this->VariableNames[60] = "SLIP_U";
  this->VariableNames[61] = "SLIP_V";
  this->VariableNames[62] = "SLIP_W";
  this->VariableNames[63] = "SDR";
  this->VariableNames[64] = "SDR_M1";
  this->VariableNames[65] = "SDR_M2";
  this->VariableNames[66] = "POLLUT4";
  this->VariableNames[67] = "GRANULAR_TEMPERATURE";
  this->VariableNames[68] = "GRANULAR_TEMPERATURE_M1";
  this->VariableNames[69] = "GRANULAR_TEMPERATURE_M2";
  this->VariableNames[70] = "VFLUX";
  this->VariableNames[71] = "V71";
  this->VariableNames[72] = "V72";
  this->VariableNames[73] = "V73";
  this->VariableNames[74] = "V74";
  this->VariableNames[75] = "V75";
  this->VariableNames[76] = "V76";
  this->VariableNames[77] = "V77";
  this->VariableNames[78] = "V78";
  this->VariableNames[79] = "V79";
  this->VariableNames[80] = "VFLUX_M1";
  this->VariableNames[81] = "V81";
  this->VariableNames[82] = "V82";
  this->VariableNames[83] = "V83";
  this->VariableNames[84] = "V84";
  this->VariableNames[85] = "V85";
  this->VariableNames[86] = "V86";
  this->VariableNames[87] = "V87";
  this->VariableNames[88] = "V88";
  this->VariableNames[89] = "V89";
  this->VariableNames[90] = "VFLUX_M2";
  this->VariableNames[91] = "DO_QNET";
  this->VariableNames[92] = "DO_QTRANS";
  this->VariableNames[93] = "DO_QREFL";
  this->VariableNames[94] = "DO_QABS";
  this->VariableNames[95] = "POLLUT5";
  this->VariableNames[96] = "WALL_DIST";
  this->VariableNames[97] = "SOLAR_SOURCE";
  this->VariableNames[98] = "SOLAR_QREFL";
  this->VariableNames[99] = "SOLAR_QABS";
  this->VariableNames[100] = "SOLAR_QTRANS";
  this->VariableNames[101] = "DENSITY";
  this->VariableNames[102] = "MU_LAM";
  this->VariableNames[103] = "MU_TURB";
  this->VariableNames[104] = "CP";
  this->VariableNames[105] = "KTC";
  this->VariableNames[106] = "VGS_DTRM";
  this->VariableNames[107] = "VGF_DTRM";
  this->VariableNames[108] = "RSTRESS";
  this->VariableNames[109] = "THREAD_RAD_FLUX";
  this->VariableNames[110] = "SPE_Q";
  this->VariableNames[111] = "X_VELOCITY";
  this->VariableNames[112] = "Y_VELOCITY";
  this->VariableNames[113] = "Z_VELOCITY";
  this->VariableNames[114] = "WALL_VELOCITY";
  this->VariableNames[115] = "X_VELOCITY_M1";
  this->VariableNames[116] = "Y_VELOCITY_M1";
  this->VariableNames[117] = "Z_VELOCITY_M1";
  this->VariableNames[118] = "PHASE_MASS";
  this->VariableNames[119] = "TKE_M1";
  this->VariableNames[120] = "TED_M1";
  this->VariableNames[121] = "POLLUT6";
  this->VariableNames[122] = "X_VELOCITY_M2";
  this->VariableNames[123] = "Y_VELOCITY_M2";
  this->VariableNames[124] = "Z_VELOCITY_M2";
  this->VariableNames[125] = "V125";
  this->VariableNames[126] = "TKE_M2";
  this->VariableNames[127] = "TED_M2";
  this->VariableNames[128] = "RUU";
  this->VariableNames[129] = "RVV";
  this->VariableNames[130] = "RWW";
  this->VariableNames[131] = "RUV";
  this->VariableNames[132] = "RVW";
  this->VariableNames[133] = "RUW";
  this->VariableNames[134] = "DPMS_EROSION";
  this->VariableNames[135] = "DPMS_ACCRETION";
  this->VariableNames[136] = "FMEAN2";
  this->VariableNames[137] = "FVAR2";
  this->VariableNames[138] = "ENTHALPY_M1";
  this->VariableNames[139] = "ENTHALPY_M2";
  this->VariableNames[140] = "FMEAN_M1";
  this->VariableNames[141] = "FMEAN_M2";
  this->VariableNames[142] = "FVAR_M1";
  this->VariableNames[143] = "FVAR_M2";
  this->VariableNames[144] = "FMEAN2_M1";
  this->VariableNames[145] = "FMEAN2_M2";
  this->VariableNames[146] = "FVAR2_M1";
  this->VariableNames[147] = "FVAR2_M2";
  this->VariableNames[148] = "PREMIXC_M1";
  this->VariableNames[149] = "PREMIXC_M2";
  this->VariableNames[150] = "VOF";
  this->VariableNames[151] = "VOF_1";
  this->VariableNames[152] = "VOF_2";
  this->VariableNames[153] = "VOF_3";
  this->VariableNames[154] = "VOF_4";
  this->VariableNames[155] = "V155";
  this->VariableNames[156] = "V156";
  this->VariableNames[157] = "V157";
  this->VariableNames[158] = "V158";
  this->VariableNames[159] = "V159";
  this->VariableNames[160] = "VOF_M1";
  this->VariableNames[161] = "VOF_1_M1";
  this->VariableNames[162] = "VOF_2_M1";
  this->VariableNames[163] = "VOF_3_M1";
  this->VariableNames[164] = "VOF_4_M1";
  this->VariableNames[165] = "V165";
  this->VariableNames[166] = "V166";
  this->VariableNames[167] = "V167";
  this->VariableNames[168] = "V168";
  this->VariableNames[169] = "V169";
  this->VariableNames[170] = "VOF_M2";
  this->VariableNames[171] = "VOF_1_M2";
  this->VariableNames[172] = "VOF_2_M2";
  this->VariableNames[173] = "VOF_3_M2";
  this->VariableNames[174] = "VOF_4_M2";
  this->VariableNames[175] = "V175";
  this->VariableNames[176] = "V176";
  this->VariableNames[177] = "V177";
  this->VariableNames[178] = "V178";
  this->VariableNames[179] = "V179";
  this->VariableNames[180] = "VOLUME_M2";
  this->VariableNames[181] = "WALL_GRID_VELOCITY";
  this->VariableNames[182] = "V182";
  this->VariableNames[183] = "V183";
  this->VariableNames[184] = "V184";
  this->VariableNames[185] = "V185";
  this->VariableNames[186] = "V186";
  this->VariableNames[187] = "V187";
  this->VariableNames[188] = "V188";
  this->VariableNames[189] = "V189";
  this->VariableNames[190] = "SV_T_AUX";
  this->VariableNames[191] = "SV_T_AP_AUX";
  this->VariableNames[192] = "TOTAL_PRESSURE";
  this->VariableNames[193] = "TOTAL_TEMPERATURE";
  this->VariableNames[194] = "NRBC_DC";
  this->VariableNames[195] = "DP_TMFR";
  this->VariableNames[196] = "V196";
  this->VariableNames[197] = "V197";
  this->VariableNames[198] = "V198";
  this->VariableNames[199] = "V199";
  this->VariableNames[200] = "SV_Y_0";
  this->VariableNames[201] = "SV_Y_1";
  this->VariableNames[202] = "SV_Y_2";
  this->VariableNames[203] = "SV_Y_3";
  this->VariableNames[204] = "SV_Y_4";
  this->VariableNames[205] = "SV_Y_5";
  this->VariableNames[206] = "SV_Y_6";
  this->VariableNames[207] = "SV_Y_7";
  this->VariableNames[208] = "SV_Y_8";
  this->VariableNames[209] = "SV_Y_9";
  this->VariableNames[210] = "SV_Y_10";
  this->VariableNames[211] = "SV_Y_11";
  this->VariableNames[212] = "SV_Y_12";
  this->VariableNames[213] = "SV_Y_13";
  this->VariableNames[214] = "SV_Y_14";
  this->VariableNames[215] = "SV_Y_15";
  this->VariableNames[216] = "SV_Y_16";
  this->VariableNames[217] = "SV_Y_17";
  this->VariableNames[218] = "SV_Y_18";
  this->VariableNames[219] = "SV_Y_19";
  this->VariableNames[220] = "SV_Y_20";
  this->VariableNames[221] = "SV_Y_21";
  this->VariableNames[222] = "SV_Y_22";
  this->VariableNames[223] = "SV_Y_23";
  this->VariableNames[224] = "SV_Y_24";
  this->VariableNames[225] = "SV_Y_25";
  this->VariableNames[226] = "SV_Y_26";
  this->VariableNames[227] = "SV_Y_27";
  this->VariableNames[228] = "SV_Y_28";
  this->VariableNames[229] = "SV_Y_29";
  this->VariableNames[230] = "SV_Y_30";
  this->VariableNames[231] = "SV_Y_31";
  this->VariableNames[232] = "SV_Y_32";
  this->VariableNames[233] = "SV_Y_33";
  this->VariableNames[234] = "SV_Y_34";
  this->VariableNames[235] = "SV_Y_35";
  this->VariableNames[236] = "SV_Y_36";
  this->VariableNames[237] = "SV_Y_37";
  this->VariableNames[238] = "SV_Y_38";
  this->VariableNames[239] = "SV_Y_39";
  this->VariableNames[240] = "SV_Y_40";
  this->VariableNames[241] = "SV_Y_41";
  this->VariableNames[242] = "SV_Y_42";
  this->VariableNames[243] = "SV_Y_43";
  this->VariableNames[244] = "SV_Y_44";
  this->VariableNames[245] = "SV_Y_45";
  this->VariableNames[246] = "SV_Y_46";
  this->VariableNames[247] = "SV_Y_47";
  this->VariableNames[248] = "SV_Y_48";
  this->VariableNames[249] = "SV_Y_49";
  this->VariableNames[250] = "SV_M1_Y_0";
  this->VariableNames[251] = "SV_M1_Y_1";
  this->VariableNames[252] = "SV_M1_Y_2";
  this->VariableNames[253] = "SV_M1_Y_3";
  this->VariableNames[254] = "SV_M1_Y_4";
  this->VariableNames[255] = "SV_M1_Y_5";
  this->VariableNames[256] = "SV_M1_Y_6";
  this->VariableNames[257] = "SV_M1_Y_7";
  this->VariableNames[258] = "SV_M1_Y_8";
  this->VariableNames[259] = "SV_M1_Y_9";
  this->VariableNames[260] = "SV_M1_Y_10";
  this->VariableNames[261] = "SV_M1_Y_11";
  this->VariableNames[262] = "SV_M1_Y_12";
  this->VariableNames[263] = "SV_M1_Y_13";
  this->VariableNames[264] = "SV_M1_Y_14";
  this->VariableNames[265] = "SV_M1_Y_15";
  this->VariableNames[266] = "SV_M1_Y_16";
  this->VariableNames[267] = "SV_M1_Y_17";
  this->VariableNames[268] = "SV_M1_Y_18";
  this->VariableNames[269] = "SV_M1_Y_19";
  this->VariableNames[270] = "SV_M1_Y_20";
  this->VariableNames[271] = "SV_M1_Y_21";
  this->VariableNames[272] = "SV_M1_Y_22";
  this->VariableNames[273] = "SV_M1_Y_23";
  this->VariableNames[274] = "SV_M1_Y_24";
  this->VariableNames[275] = "SV_M1_Y_25";
  this->VariableNames[276] = "SV_M1_Y_26";
  this->VariableNames[277] = "SV_M1_Y_27";
  this->VariableNames[278] = "SV_M1_Y_28";
  this->VariableNames[279] = "SV_M1_Y_29";
  this->VariableNames[280] = "SV_M1_Y_30";
  this->VariableNames[281] = "SV_M1_Y_31";
  this->VariableNames[282] = "SV_M1_Y_32";
  this->VariableNames[283] = "SV_M1_Y_33";
  this->VariableNames[284] = "SV_M1_Y_34";
  this->VariableNames[285] = "SV_M1_Y_35";
  this->VariableNames[286] = "SV_M1_Y_36";
  this->VariableNames[287] = "SV_M1_Y_37";
  this->VariableNames[288] = "SV_M1_Y_38";
  this->VariableNames[289] = "SV_M1_Y_39";
  this->VariableNames[290] = "SV_M1_Y_40";
  this->VariableNames[291] = "SV_M1_Y_41";
  this->VariableNames[292] = "SV_M1_Y_42";
  this->VariableNames[293] = "SV_M1_Y_43";
  this->VariableNames[294] = "SV_M1_Y_44";
  this->VariableNames[295] = "SV_M1_Y_45";
  this->VariableNames[296] = "SV_M1_Y_46";
  this->VariableNames[297] = "SV_M1_Y_47";
  this->VariableNames[298] = "SV_M1_Y_48";
  this->VariableNames[299] = "SV_M1_Y_49";
  this->VariableNames[300] = "SV_M2_Y_0";
  this->VariableNames[301] = "SV_M2_Y_1";
  this->VariableNames[302] = "SV_M2_Y_2";
  this->VariableNames[303] = "SV_M2_Y_3";
  this->VariableNames[304] = "SV_M2_Y_4";
  this->VariableNames[305] = "SV_M2_Y_5";
  this->VariableNames[306] = "SV_M2_Y_6";
  this->VariableNames[307] = "SV_M2_Y_7";
  this->VariableNames[308] = "SV_M2_Y_8";
  this->VariableNames[309] = "SV_M2_Y_9";
  this->VariableNames[310] = "SV_M2_Y_10";
  this->VariableNames[311] = "SV_M2_Y_11";
  this->VariableNames[312] = "SV_M2_Y_12";
  this->VariableNames[313] = "SV_M2_Y_13";
  this->VariableNames[314] = "SV_M2_Y_14";
  this->VariableNames[315] = "SV_M2_Y_15";
  this->VariableNames[316] = "SV_M2_Y_16";
  this->VariableNames[317] = "SV_M2_Y_17";
  this->VariableNames[318] = "SV_M2_Y_18";
  this->VariableNames[319] = "SV_M2_Y_19";
  this->VariableNames[320] = "SV_M2_Y_20";
  this->VariableNames[321] = "SV_M2_Y_21";
  this->VariableNames[322] = "SV_M2_Y_22";
  this->VariableNames[323] = "SV_M2_Y_23";
  this->VariableNames[324] = "SV_M2_Y_24";
  this->VariableNames[325] = "SV_M2_Y_25";
  this->VariableNames[326] = "SV_M2_Y_26";
  this->VariableNames[327] = "SV_M2_Y_27";
  this->VariableNames[328] = "SV_M2_Y_28";
  this->VariableNames[329] = "SV_M2_Y_29";
  this->VariableNames[330] = "SV_M2_Y_30";
  this->VariableNames[331] = "SV_M2_Y_31";
  this->VariableNames[332] = "SV_M2_Y_32";
  this->VariableNames[333] = "SV_M2_Y_33";
  this->VariableNames[334] = "SV_M2_Y_34";
  this->VariableNames[335] = "SV_M2_Y_35";
  this->VariableNames[336] = "SV_M2_Y_36";
  this->VariableNames[337] = "SV_M2_Y_37";
  this->VariableNames[338] = "SV_M2_Y_38";
  this->VariableNames[339] = "SV_M2_Y_39";
  this->VariableNames[340] = "SV_M2_Y_40";
  this->VariableNames[341] = "SV_M2_Y_41";
  this->VariableNames[342] = "SV_M2_Y_42";
  this->VariableNames[343] = "SV_M2_Y_43";
  this->VariableNames[344] = "SV_M2_Y_44";
  this->VariableNames[345] = "SV_M2_Y_45";
  this->VariableNames[346] = "SV_M2_Y_46";
  this->VariableNames[347] = "SV_M2_Y_47";
  this->VariableNames[348] = "SV_M2_Y_48";
  this->VariableNames[349] = "SV_M2_Y_49";
  this->VariableNames[350] = "DR_SURF_0";
  this->VariableNames[351] = "DR_SURF_1";
  this->VariableNames[352] = "DR_SURF_2";
  this->VariableNames[353] = "DR_SURF_3";
  this->VariableNames[354] = "DR_SURF_4";
  this->VariableNames[355] = "DR_SURF_5";
  this->VariableNames[356] = "DR_SURF_6";
  this->VariableNames[357] = "DR_SURF_7";
  this->VariableNames[358] = "DR_SURF_8";
  this->VariableNames[359] = "DR_SURF_9";
  this->VariableNames[360] = "DR_SURF_10";
  this->VariableNames[361] = "DR_SURF_11";
  this->VariableNames[362] = "DR_SURF_12";
  this->VariableNames[363] = "DR_SURF_13";
  this->VariableNames[364] = "DR_SURF_14";
  this->VariableNames[365] = "DR_SURF_15";
  this->VariableNames[366] = "DR_SURF_16";
  this->VariableNames[367] = "DR_SURF_17";
  this->VariableNames[368] = "DR_SURF_18";
  this->VariableNames[369] = "DR_SURF_19";
  this->VariableNames[370] = "DR_SURF_20";
  this->VariableNames[371] = "DR_SURF_21";
  this->VariableNames[372] = "DR_SURF_22";
  this->VariableNames[373] = "DR_SURF_23";
  this->VariableNames[374] = "DR_SURF_24";
  this->VariableNames[375] = "DR_SURF_25";
  this->VariableNames[376] = "DR_SURF_26";
  this->VariableNames[377] = "DR_SURF_27";
  this->VariableNames[378] = "DR_SURF_28";
  this->VariableNames[379] = "DR_SURF_29";
  this->VariableNames[380] = "DR_SURF_30";
  this->VariableNames[381] = "DR_SURF_31";
  this->VariableNames[382] = "DR_SURF_32";
  this->VariableNames[383] = "DR_SURF_33";
  this->VariableNames[384] = "DR_SURF_34";
  this->VariableNames[385] = "DR_SURF_35";
  this->VariableNames[386] = "DR_SURF_36";
  this->VariableNames[387] = "DR_SURF_37";
  this->VariableNames[388] = "DR_SURF_38";
  this->VariableNames[389] = "DR_SURF_39";
  this->VariableNames[390] = "DR_SURF_40";
  this->VariableNames[391] = "DR_SURF_41";
  this->VariableNames[392] = "DR_SURF_42";
  this->VariableNames[393] = "DR_SURF_43";
  this->VariableNames[394] = "DR_SURF_44";
  this->VariableNames[395] = "DR_SURF_45";
  this->VariableNames[396] = "DR_SURF_46";
  this->VariableNames[397] = "DR_SURF_47";
  this->VariableNames[398] = "DR_SURF_48";
  this->VariableNames[399] = "DR_SURF_49";
  this->VariableNames[400] = "PRESSURE_MEAN";
  this->VariableNames[401] = "PRESSURE_RMS";
  this->VariableNames[402] = "X_VELOCITY_MEAN";
  this->VariableNames[403] = "X_VELOCITY_RMS";
  this->VariableNames[404] = "Y_VELOCITY_MEAN";
  this->VariableNames[405] = "Y_VELOCITY_RMS";
  this->VariableNames[406] = "Z_VELOCITY_MEAN";
  this->VariableNames[407] = "Z_VELOCITY_RMS";
  this->VariableNames[408] = "TEMPERATURE_MEAN";
  this->VariableNames[409] = "TEMPERATURE_RMS";
  this->VariableNames[410] = "VOF_MEAN";
  this->VariableNames[411] = "VOF_RMS";
  this->VariableNames[412] = "PRESSURE_M1";
  this->VariableNames[413] = "PRESSURE_M2";
  this->VariableNames[414] = "GRANULAR_TEMPERATURE_MEAN";
  this->VariableNames[415] = "GRANULAR_TEMPERATURE_RMS";
  this->VariableNames[416] = "V416";
  this->VariableNames[417] = "V417";
  this->VariableNames[418] = "V418";
  this->VariableNames[419] = "V419";
  this->VariableNames[420] = "V420";
  this->VariableNames[421] = "V421";
  this->VariableNames[422] = "V422";
  this->VariableNames[423] = "V423";
  this->VariableNames[424] = "V424";
  this->VariableNames[425] = "V425";
  this->VariableNames[426] = "V426";
  this->VariableNames[427] = "V427";
  this->VariableNames[428] = "V428";
  this->VariableNames[429] = "V429";
  this->VariableNames[430] = "V430";
  this->VariableNames[431] = "V431";
  this->VariableNames[432] = "V432";
  this->VariableNames[433] = "V433";
  this->VariableNames[434] = "V434";
  this->VariableNames[435] = "V435";
  this->VariableNames[436] = "V436";
  this->VariableNames[437] = "V437";
  this->VariableNames[438] = "V438";
  this->VariableNames[439] = "V439";
  this->VariableNames[440] = "V440";
  this->VariableNames[441] = "V441";
  this->VariableNames[442] = "V442";
  this->VariableNames[443] = "V443";
  this->VariableNames[444] = "V444";
  this->VariableNames[445] = "V445";
  this->VariableNames[446] = "V446";
  this->VariableNames[447] = "V447";
  this->VariableNames[448] = "V448";
  this->VariableNames[449] = "V449";
  this->VariableNames[450] = "DPMS_Y_0";
  this->VariableNames[451] = "DPMS_Y_1";
  this->VariableNames[452] = "DPMS_Y_2";
  this->VariableNames[453] = "DPMS_Y_3";
  this->VariableNames[454] = "DPMS_Y_4";
  this->VariableNames[455] = "DPMS_Y_5";
  this->VariableNames[456] = "DPMS_Y_6";
  this->VariableNames[457] = "DPMS_Y_7";
  this->VariableNames[458] = "DPMS_Y_8";
  this->VariableNames[459] = "DPMS_Y_9";
  this->VariableNames[460] = "DPMS_Y_10";
  this->VariableNames[461] = "DPMS_Y_11";
  this->VariableNames[462] = "DPMS_Y_12";
  this->VariableNames[463] = "DPMS_Y_13";
  this->VariableNames[464] = "DPMS_Y_14";
  this->VariableNames[465] = "DPMS_Y_15";
  this->VariableNames[466] = "DPMS_Y_16";
  this->VariableNames[467] = "DPMS_Y_17";
  this->VariableNames[468] = "DPMS_Y_18";
  this->VariableNames[469] = "DPMS_Y_19";
  this->VariableNames[470] = "DPMS_Y_20";
  this->VariableNames[471] = "DPMS_Y_21";
  this->VariableNames[472] = "DPMS_Y_22";
  this->VariableNames[473] = "DPMS_Y_23";
  this->VariableNames[474] = "DPMS_Y_24";
  this->VariableNames[475] = "DPMS_Y_25";
  this->VariableNames[476] = "DPMS_Y_26";
  this->VariableNames[477] = "DPMS_Y_27";
  this->VariableNames[478] = "DPMS_Y_28";
  this->VariableNames[479] = "DPMS_Y_29";
  this->VariableNames[480] = "DPMS_Y_30";
  this->VariableNames[481] = "DPMS_Y_31";
  this->VariableNames[482] = "DPMS_Y_32";
  this->VariableNames[483] = "DPMS_Y_33";
  this->VariableNames[484] = "DPMS_Y_34";
  this->VariableNames[485] = "DPMS_Y_35";
  this->VariableNames[486] = "DPMS_Y_36";
  this->VariableNames[487] = "DPMS_Y_37";
  this->VariableNames[488] = "DPMS_Y_38";
  this->VariableNames[489] = "DPMS_Y_39";
  this->VariableNames[490] = "DPMS_Y_40";
  this->VariableNames[491] = "DPMS_Y_41";
  this->VariableNames[492] = "DPMS_Y_42";
  this->VariableNames[493] = "DPMS_Y_43";
  this->VariableNames[494] = "DPMS_Y_44";
  this->VariableNames[495] = "DPMS_Y_45";
  this->VariableNames[496] = "DPMS_Y_46";
  this->VariableNames[497] = "DPMS_Y_47";
  this->VariableNames[498] = "DPMS_Y_48";
  this->VariableNames[499] = "DPMS_Y_49";
  this->VariableNames[500] = "NUT";
  this->VariableNames[501] = "NUT_M1";
  this->VariableNames[502] = "NUT_M2";
  this->VariableNames[503] = "RUU_M1";
  this->VariableNames[504] = "RVV_M1";
  this->VariableNames[505] = "RWW_M1";
  this->VariableNames[506] = "RUV_M1";
  this->VariableNames[507] = "RVW_M1";
  this->VariableNames[508] = "RUW_M1";
  this->VariableNames[509] = "RUU_M2";
  this->VariableNames[510] = "RVV_M2";
  this->VariableNames[511] = "RWW_M2";
  this->VariableNames[512] = "RUV_M2";
  this->VariableNames[513] = "RVW_M2";
  this->VariableNames[514] = "RUW_M2";
  this->VariableNames[515] = "ENERGY_M1";
  this->VariableNames[516] = "ENERGY_M2";
  this->VariableNames[517] = "DENSITY_M1";
  this->VariableNames[518] = "DENSITY_M2";
  this->VariableNames[519] = "DPMS_PDF_1";
  this->VariableNames[520] = "DPMS_PDF_2";
  this->VariableNames[521] = "V2";
  this->VariableNames[522] = "V2_M1";
  this->VariableNames[523] = "V2_M2";
  this->VariableNames[524] = "FEL";
  this->VariableNames[525] = "FEL_M1";
  this->VariableNames[526] = "FEL_M2";
  this->VariableNames[527] = "V527";
  this->VariableNames[528] = "V528";
  this->VariableNames[529] = "V529";
  this->VariableNames[530] = "SHELL_CELL_T";
  this->VariableNames[531] = "SHELL_FACE_T";
  this->VariableNames[532] = "V532";
  this->VariableNames[533] = "V533";
  this->VariableNames[534] = "V534";
  this->VariableNames[535] = "V535";
  this->VariableNames[536] = "V536";
  this->VariableNames[537] = "V537";
  this->VariableNames[538] = "V538";
  this->VariableNames[539] = "V539";
  this->VariableNames[540] = "DPMS_TKE";
  this->VariableNames[541] = "DPMS_D";
  this->VariableNames[542] = "DPMS_O";
  this->VariableNames[543] = "DPMS_TKE_RUU";
  this->VariableNames[544] = "DPMS_TKE_RVV";
  this->VariableNames[545] = "DPMS_TKE_RWW";
  this->VariableNames[546] = "DPMS_TKE_RUV";
  this->VariableNames[547] = "DPMS_TKE_RVW";
  this->VariableNames[548] = "DPMS_TKE_RUW";
  this->VariableNames[549] = "DPMS_DS_MASS";
  this->VariableNames[550] = "DPMS_DS_ENERGY";
  this->VariableNames[551] = "DPMS_DS_TKE";
  this->VariableNames[552] = "DPMS_DS_D";
  this->VariableNames[553] = "DPMS_DS_O";
  this->VariableNames[554] = "DPMS_DS_TKE_RUU";
  this->VariableNames[555] = "DPMS_DS_TKE_RVV";
  this->VariableNames[556] = "DPMS_DS_TKE_RWW";
  this->VariableNames[557] = "DPMS_DS_TKE_RUV";
  this->VariableNames[558] = "DPMS_DS_TKE_RVW";
  this->VariableNames[559] = "DPMS_DS_TKE_RUW";
  this->VariableNames[560] = "DPMS_DS_PDF_1";
  this->VariableNames[561] = "DPMS_DS_PDF_2";
  this->VariableNames[562] = "DPMS_DS_EMISS";
  this->VariableNames[563] = "DPMS_DS_ABS";
  this->VariableNames[564] = "DPMS_DS_SCAT";
  this->VariableNames[565] = "DPMS_DS_BURNOUT";
  this->VariableNames[566] = "DPMS_DS_MOM";
  this->VariableNames[567] = "DPMS_DS_WSWIRL";
  this->VariableNames[568] = "V568";
  this->VariableNames[569] = "V569";
  this->VariableNames[570] = "V570";
  this->VariableNames[571] = "V571";
  this->VariableNames[572] = "V572";
  this->VariableNames[573] = "V573";
  this->VariableNames[574] = "V574";
  this->VariableNames[575] = "V575";
  this->VariableNames[576] = "V576";
  this->VariableNames[577] = "V577";
  this->VariableNames[578] = "V578";
  this->VariableNames[579] = "V579";
  this->VariableNames[580] = "V570";
  this->VariableNames[581] = "V581";
  this->VariableNames[582] = "V582";
  this->VariableNames[583] = "V583";
  this->VariableNames[584] = "V584";
  this->VariableNames[585] = "V585";
  this->VariableNames[586] = "V586";
  this->VariableNames[587] = "V587";
  this->VariableNames[588] = "V588";
  this->VariableNames[589] = "V589";
  this->VariableNames[590] = "V590";
  this->VariableNames[591] = "V591";
  this->VariableNames[592] = "V592";
  this->VariableNames[593] = "V593";
  this->VariableNames[594] = "V594";
  this->VariableNames[595] = "V595";
  this->VariableNames[596] = "V596";
  this->VariableNames[597] = "V597";
  this->VariableNames[598] = "V598";
  this->VariableNames[599] = "V599";
  this->VariableNames[600] = "DELH";
  this->VariableNames[601] = "DPMS_MOM_AP";
  this->VariableNames[602] = "DPMS_WSWIRL_AP";
  this->VariableNames[603] = "X_PULL";
  this->VariableNames[604] = "Y_PULL";
  this->VariableNames[605] = "Z_PULL";
  this->VariableNames[606] = "LIQF";
  this->VariableNames[607] = "V607";
  this->VariableNames[608] = "V608";
  this->VariableNames[609] = "V609";
  this->VariableNames[610] = "PDFT_QBAR";
  this->VariableNames[611] = "PDFT_PHI";
  this->VariableNames[612] = "PDFT_Q_TA";
  this->VariableNames[613] = "PDFT_SVOL_TA";
  this->VariableNames[614] = "PDFT_MASS_TA";
  this->VariableNames[615] = "PDFT_T4_TA";
  this->VariableNames[616] = "V616";
  this->VariableNames[617] = "V617";
  this->VariableNames[618] = "V618";
  this->VariableNames[619] = "V619";
  this->VariableNames[620] = "V620";
  this->VariableNames[621] = "V621";
  this->VariableNames[622] = "V622";
  this->VariableNames[623] = "V623";
  this->VariableNames[624] = "V624";
  this->VariableNames[625] = "V625";
  this->VariableNames[626] = "V626";
  this->VariableNames[627] = "V627";
  this->VariableNames[628] = "V628";
  this->VariableNames[629] = "V629";
  this->VariableNames[630] = "SCAD_LES";
  this->VariableNames[631] = "V631";
  this->VariableNames[632] = "V632";
  this->VariableNames[633] = "V633";
  this->VariableNames[634] = "V634";
  this->VariableNames[635] = "V635";
  this->VariableNames[636] = "V636";
  this->VariableNames[637] = "V637";
  this->VariableNames[638] = "V638";
  this->VariableNames[639] = "V639";
  this->VariableNames[640] = "V640";
  this->VariableNames[641] = "V641";
  this->VariableNames[642] = "V642";
  this->VariableNames[643] = "V643";
  this->VariableNames[644] = "V644";
  this->VariableNames[645] = "CREV_MASS";
  this->VariableNames[646] = "CREV_ENRG";
  this->VariableNames[647] = "CREV_MOM";
  this->VariableNames[648] = "V648";
  this->VariableNames[649] = "V649";
  this->VariableNames[650] = "XF_ACOUSTICS_MODEL";
  this->VariableNames[651] = "XF_RF_AC_RECEIVERS_DATA";
  this->VariableNames[652] = "SV_DPDT_RMS";
  this->VariableNames[653] = "SV_PRESSURE_M1";
  this->VariableNames[654] = "XF_RF_AC_PERIODIC_INDEX";
  this->VariableNames[655] = "XF_RF_AC_PERIODIC_PS";
  this->VariableNames[656] = "XF_RF_AC_F_NORMAL";
  this->VariableNames[657] = "XF_RF_AC_F_CENTROID";
  this->VariableNames[658] = "V658";
  this->VariableNames[659] = "V659";
  this->VariableNames[660] = "IGNITE";
  this->VariableNames[661] = "IGNITE_M1";
  this->VariableNames[662] = "IGNITE_M2";
  this->VariableNames[663] = "IGNITE_RATE";
  this->VariableNames[664] = "V664";
  this->VariableNames[665] = "V665";
  this->VariableNames[666] = "V666";
  this->VariableNames[667] = "V667";
  this->VariableNames[668] = "V668";
  this->VariableNames[669] = "V669";
  this->VariableNames[670] = "V670";
  this->VariableNames[671] = "V671";
  this->VariableNames[672] = "V672";
  this->VariableNames[673] = "V673";
  this->VariableNames[674] = "V674";
  this->VariableNames[675] = "V675";
  this->VariableNames[676] = "V676";
  this->VariableNames[677] = "V677";
  this->VariableNames[678] = "V678";
  this->VariableNames[679] = "V679";
  this->VariableNames[680] = "V680";
  this->VariableNames[681] = "V681";
  this->VariableNames[682] = "V682";
  this->VariableNames[683] = "V683";
  this->VariableNames[684] = "V684";
  this->VariableNames[685] = "V685";
  this->VariableNames[686] = "V686";
  this->VariableNames[687] = "V687";
  this->VariableNames[688] = "V688";
  this->VariableNames[689] = "V689";
  this->VariableNames[690] = "V690";
  this->VariableNames[691] = "V691";
  this->VariableNames[692] = "V692";
  this->VariableNames[693] = "V693";
  this->VariableNames[694] = "V694";
  this->VariableNames[695] = "V695";
  this->VariableNames[696] = "V696";
  this->VariableNames[697] = "V697";
  this->VariableNames[698] = "V698";
  this->VariableNames[699] = "V699";
  this->VariableNames[700] = "UDS_0";
  this->VariableNames[701] = "UDS_1";
  this->VariableNames[702] = "UDS_2";
  this->VariableNames[703] = "UDS_3";
  this->VariableNames[704] = "UDS_4";
  this->VariableNames[705] = "UDS_5";
  this->VariableNames[706] = "UDS_6";
  this->VariableNames[707] = "UDS_7";
  this->VariableNames[708] = "UDS_8";
  this->VariableNames[709] = "UDS_9";
  this->VariableNames[710] = "UDS_10";
  this->VariableNames[711] = "UDS_11";
  this->VariableNames[712] = "UDS_12";
  this->VariableNames[713] = "UDS_13";
  this->VariableNames[714] = "UDS_14";
  this->VariableNames[715] = "UDS_15";
  this->VariableNames[716] = "UDS_16";
  this->VariableNames[717] = "UDS_17";
  this->VariableNames[718] = "UDS_18";
  this->VariableNames[719] = "UDS_19";
  this->VariableNames[720] = "UDS_20";
  this->VariableNames[721] = "UDS_21";
  this->VariableNames[722] = "UDS_22";
  this->VariableNames[723] = "UDS_23";
  this->VariableNames[724] = "UDS_24";
  this->VariableNames[725] = "UDS_25";
  this->VariableNames[726] = "UDS_26";
  this->VariableNames[727] = "UDS_27";
  this->VariableNames[728] = "UDS_28";
  this->VariableNames[729] = "UDS_29";
  this->VariableNames[730] = "UDS_30";
  this->VariableNames[731] = "UDS_31";
  this->VariableNames[732] = "UDS_32";
  this->VariableNames[733] = "UDS_33";
  this->VariableNames[734] = "UDS_34";
  this->VariableNames[735] = "UDS_35";
  this->VariableNames[736] = "UDS_36";
  this->VariableNames[737] = "UDS_37";
  this->VariableNames[738] = "UDS_38";
  this->VariableNames[739] = "UDS_39";
  this->VariableNames[740] = "UDS_40";
  this->VariableNames[741] = "UDS_41";
  this->VariableNames[742] = "UDS_42";
  this->VariableNames[743] = "UDS_43";
  this->VariableNames[744] = "UDS_44";
  this->VariableNames[745] = "UDS_45";
  this->VariableNames[746] = "UDS_46";
  this->VariableNames[747] = "UDS_47";
  this->VariableNames[748] = "UDS_48";
  this->VariableNames[749] = "UDS_49";
  this->VariableNames[750] = "UDS_M1_0";
  this->VariableNames[751] = "UDS_M1_1";
  this->VariableNames[752] = "UDS_M1_2";
  this->VariableNames[753] = "UDS_M1_3";
  this->VariableNames[754] = "UDS_M1_4";
  this->VariableNames[755] = "UDS_M1_5";
  this->VariableNames[756] = "UDS_M1_6";
  this->VariableNames[757] = "UDS_M1_7";
  this->VariableNames[758] = "UDS_M1_8";
  this->VariableNames[759] = "UDS_M1_9";
  this->VariableNames[760] = "UDS_M1_10";
  this->VariableNames[761] = "UDS_M1_11";
  this->VariableNames[762] = "UDS_M1_12";
  this->VariableNames[763] = "UDS_M1_13";
  this->VariableNames[764] = "UDS_M1_14";
  this->VariableNames[765] = "UDS_M1_15";
  this->VariableNames[766] = "UDS_M1_16";
  this->VariableNames[767] = "UDS_M1_17";
  this->VariableNames[768] = "UDS_M1_18";
  this->VariableNames[769] = "UDS_M1_19";
  this->VariableNames[770] = "UDS_M1_20";
  this->VariableNames[771] = "UDS_M1_21";
  this->VariableNames[772] = "UDS_M1_22";
  this->VariableNames[773] = "UDS_M1_23";
  this->VariableNames[774] = "UDS_M1_24";
  this->VariableNames[775] = "UDS_M1_25";
  this->VariableNames[776] = "UDS_M1_26";
  this->VariableNames[777] = "UDS_M1_27";
  this->VariableNames[778] = "UDS_M1_28";
  this->VariableNames[779] = "UDS_M1_29";
  this->VariableNames[780] = "UDS_M1_30";
  this->VariableNames[781] = "UDS_M1_31";
  this->VariableNames[782] = "UDS_M1_32";
  this->VariableNames[783] = "UDS_M1_33";
  this->VariableNames[784] = "UDS_M1_34";
  this->VariableNames[785] = "UDS_M1_35";
  this->VariableNames[786] = "UDS_M1_36";
  this->VariableNames[787] = "UDS_M1_37";
  this->VariableNames[788] = "UDS_M1_38";
  this->VariableNames[789] = "UDS_M1_39";
  this->VariableNames[790] = "UDS_M1_40";
  this->VariableNames[791] = "UDS_M1_41";
  this->VariableNames[792] = "UDS_M1_42";
  this->VariableNames[793] = "UDS_M1_43";
  this->VariableNames[794] = "UDS_M1_44";
  this->VariableNames[795] = "UDS_M1_45";
  this->VariableNames[796] = "UDS_M1_46";
  this->VariableNames[797] = "UDS_M1_47";
  this->VariableNames[798] = "UDS_M1_48";
  this->VariableNames[799] = "UDS_M1_49";
  this->VariableNames[800] = "UDS_M2_0";
  this->VariableNames[801] = "UDS_M2_1";
  this->VariableNames[802] = "UDS_M2_2";
  this->VariableNames[803] = "UDS_M2_3";
  this->VariableNames[804] = "UDS_M2_4";
  this->VariableNames[805] = "UDS_M2_5";
  this->VariableNames[806] = "UDS_M2_6";
  this->VariableNames[807] = "UDS_M2_7";
  this->VariableNames[808] = "UDS_M2_8";
  this->VariableNames[809] = "UDS_M2_9";
  this->VariableNames[810] = "UDS_M2_10";
  this->VariableNames[811] = "UDS_M2_11";
  this->VariableNames[812] = "UDS_M2_12";
  this->VariableNames[813] = "UDS_M2_13";
  this->VariableNames[814] = "UDS_M2_14";
  this->VariableNames[815] = "UDS_M2_15";
  this->VariableNames[816] = "UDS_M2_16";
  this->VariableNames[817] = "UDS_M2_17";
  this->VariableNames[818] = "UDS_M2_18";
  this->VariableNames[819] = "UDS_M2_19";
  this->VariableNames[820] = "UDS_M2_20";
  this->VariableNames[821] = "UDS_M2_21";
  this->VariableNames[822] = "UDS_M2_22";
  this->VariableNames[823] = "UDS_M2_23";
  this->VariableNames[824] = "UDS_M2_24";
  this->VariableNames[825] = "UDS_M2_25";
  this->VariableNames[826] = "UDS_M2_26";
  this->VariableNames[827] = "UDS_M2_27";
  this->VariableNames[828] = "UDS_M2_28";
  this->VariableNames[829] = "UDS_M2_29";
  this->VariableNames[830] = "UDS_M2_30";
  this->VariableNames[831] = "UDS_M2_31";
  this->VariableNames[832] = "UDS_M2_32";
  this->VariableNames[833] = "UDS_M2_33";
  this->VariableNames[834] = "UDS_M2_34";
  this->VariableNames[835] = "UDS_M2_35";
  this->VariableNames[836] = "UDS_M2_36";
  this->VariableNames[837] = "UDS_M2_37";
  this->VariableNames[838] = "UDS_M2_38";
  this->VariableNames[839] = "UDS_M2_39";
  this->VariableNames[840] = "UDS_M2_40";
  this->VariableNames[841] = "UDS_M2_41";
  this->VariableNames[842] = "UDS_M2_42";
  this->VariableNames[843] = "UDS_M2_43";
  this->VariableNames[844] = "UDS_M2_44";
  this->VariableNames[845] = "UDS_M2_45";
  this->VariableNames[846] = "UDS_M2_46";
  this->VariableNames[847] = "UDS_M2_47";
  this->VariableNames[848] = "UDS_M2_48";
  this->VariableNames[849] = "UDS_M2_49";
  this->VariableNames[850] = "DPMS_DS_Y_0";
  this->VariableNames[851] = "DPMS_DS_Y_1";
  this->VariableNames[852] = "DPMS_DS_Y_2";
  this->VariableNames[853] = "DPMS_DS_Y_3";
  this->VariableNames[854] = "DPMS_DS_Y_4";
  this->VariableNames[855] = "DPMS_DS_Y_5";
  this->VariableNames[856] = "DPMS_DS_Y_6";
  this->VariableNames[857] = "DPMS_DS_Y_7";
  this->VariableNames[858] = "DPMS_DS_Y_8";
  this->VariableNames[859] = "DPMS_DS_Y_9";
  this->VariableNames[860] = "DPMS_DS_Y_10";
  this->VariableNames[861] = "DPMS_DS_Y_11";
  this->VariableNames[862] = "DPMS_DS_Y_12";
  this->VariableNames[863] = "DPMS_DS_Y_13";
  this->VariableNames[864] = "DPMS_DS_Y_14";
  this->VariableNames[865] = "DPMS_DS_Y_15";
  this->VariableNames[866] = "DPMS_DS_Y_16";
  this->VariableNames[867] = "DPMS_DS_Y_17";
  this->VariableNames[868] = "DPMS_DS_Y_18";
  this->VariableNames[869] = "DPMS_DS_Y_19";
  this->VariableNames[870] = "DPMS_DS_Y_20";
  this->VariableNames[871] = "DPMS_DS_Y_21";
  this->VariableNames[872] = "DPMS_DS_Y_22";
  this->VariableNames[873] = "DPMS_DS_Y_23";
  this->VariableNames[874] = "DPMS_DS_Y_24";
  this->VariableNames[875] = "DPMS_DS_Y_25";
  this->VariableNames[876] = "DPMS_DS_Y_26";
  this->VariableNames[877] = "DPMS_DS_Y_27";
  this->VariableNames[878] = "DPMS_DS_Y_28";
  this->VariableNames[879] = "DPMS_DS_Y_29";
  this->VariableNames[880] = "DPMS_DS_Y_30";
  this->VariableNames[881] = "DPMS_DS_Y_31";
  this->VariableNames[882] = "DPMS_DS_Y_32";
  this->VariableNames[883] = "DPMS_DS_Y_33";
  this->VariableNames[884] = "DPMS_DS_Y_34";
  this->VariableNames[885] = "DPMS_DS_Y_35";
  this->VariableNames[886] = "DPMS_DS_Y_36";
  this->VariableNames[887] = "DPMS_DS_Y_37";
  this->VariableNames[888] = "DPMS_DS_Y_38";
  this->VariableNames[889] = "DPMS_DS_Y_39";
  this->VariableNames[890] = "DPMS_DS_Y_40";
  this->VariableNames[891] = "DPMS_DS_Y_41";
  this->VariableNames[892] = "DPMS_DS_Y_42";
  this->VariableNames[893] = "DPMS_DS_Y_43";
  this->VariableNames[894] = "DPMS_DS_Y_44";
  this->VariableNames[895] = "DPMS_DS_Y_45";
  this->VariableNames[896] = "DPMS_DS_Y_46";
  this->VariableNames[897] = "DPMS_DS_Y_47";
  this->VariableNames[898] = "DPMS_DS_Y_48";
  this->VariableNames[899] = "DPMS_DS_Y_49";
  this->VariableNames[900] = "V900";
  this->VariableNames[901] = "V901";
  this->VariableNames[902] = "V902";
  this->VariableNames[903] = "V903";
  this->VariableNames[904] = "V904";
  this->VariableNames[905] = "V905";
  this->VariableNames[906] = "V906";
  this->VariableNames[907] = "V907";
  this->VariableNames[908] = "V908";
  this->VariableNames[909] = "V909";
  this->VariableNames[910] = "GRANULAR_PRESSURE";
  this->VariableNames[911] = "DPMS_DS_P1_S";
  this->VariableNames[912] = "DPMS_DS_P1_AP";
  this->VariableNames[913] = "DPMS_DS_P1_DIFF";
  this->VariableNames[914] = "V914";
  this->VariableNames[915] = "V915";
  this->VariableNames[916] = "V916";
  this->VariableNames[917] = "V917";
  this->VariableNames[918] = "V918";
  this->VariableNames[919] = "V919";
  this->VariableNames[920] = "V920";
  this->VariableNames[921] = "V921";
  this->VariableNames[922] = "V922";
  this->VariableNames[923] = "V923";
  this->VariableNames[924] = "V924";
  this->VariableNames[925] = "V925";
  this->VariableNames[926] = "V926";
  this->VariableNames[927] = "V927";
  this->VariableNames[928] = "V928";
  this->VariableNames[929] = "V929";
  this->VariableNames[930] = "V930";
  this->VariableNames[931] = "V931";
  this->VariableNames[932] = "V932";
  this->VariableNames[933] = "V933";
  this->VariableNames[934] = "V934";
  this->VariableNames[935] = "V935";
  this->VariableNames[936] = "V936";
  this->VariableNames[937] = "V937";
  this->VariableNames[938] = "V938";
  this->VariableNames[939] = "V939";
  this->VariableNames[940] = "V940";
  this->VariableNames[941] = "V941";
  this->VariableNames[942] = "V942";
  this->VariableNames[943] = "V943";
  this->VariableNames[944] = "V944";
  this->VariableNames[945] = "V945";
  this->VariableNames[946] = "V946";
  this->VariableNames[947] = "V947";
  this->VariableNames[948] = "V948";
  this->VariableNames[949] = "V949";
  this->VariableNames[950] = "V950";
  this->VariableNames[951] = "V951";
  this->VariableNames[952] = "V952";
  this->VariableNames[953] = "V953";
  this->VariableNames[954] = "V954";
  this->VariableNames[955] = "V955";
  this->VariableNames[956] = "V956";
  this->VariableNames[957] = "V957";
  this->VariableNames[958] = "V958";
  this->VariableNames[959] = "V959";
  this->VariableNames[960] = "V960";
  this->VariableNames[961] = "V961";
  this->VariableNames[962] = "V962";
  this->VariableNames[963] = "V963";
  this->VariableNames[964] = "V964";
  this->VariableNames[965] = "V965";
  this->VariableNames[966] = "V966";
  this->VariableNames[967] = "V967";
  this->VariableNames[968] = "V968";
  this->VariableNames[969] = "V969";
  this->VariableNames[970] = "UDM_I";
  this->VariableNames[971] = "V971";
  this->VariableNames[972] = "V972";
  this->VariableNames[973] = "V973";
  this->VariableNames[974] = "V974";
  this->VariableNames[975] = "V975";
  this->VariableNames[976] = "V976";
  this->VariableNames[977] = "V977";
  this->VariableNames[978] = "V978";
  this->VariableNames[979] = "V979";
  this->VariableNames[980] = "V980";
  this->VariableNames[981] = "V981";
  this->VariableNames[982] = "V982";
  this->VariableNames[983] = "V983";
  this->VariableNames[984] = "V984";
  this->VariableNames[985] = "V985";
  this->VariableNames[986] = "V986";
  this->VariableNames[987] = "V987";
  this->VariableNames[988] = "V988";
  this->VariableNames[989] = "V989";
  this->VariableNames[990] = "V990";
  this->VariableNames[991] = "V991";
  this->VariableNames[992] = "V992";
  this->VariableNames[993] = "V993";
  this->VariableNames[994] = "V994";
  this->VariableNames[995] = "V995";
  this->VariableNames[996] = "V996";
  this->VariableNames[997] = "V997";
  this->VariableNames[998] = "V998";
  this->VariableNames[999] = "V999";
  this->VariableNames[1000] = "V1000";
  this->VariableNames[1001] = "V1001";
  this->VariableNames[1002] = "V1002";
  this->VariableNames[1003] = "V1003";
  this->VariableNames[1004] = "V1004";
  this->VariableNames[1005] = "V1005";
  this->VariableNames[1006] = "V1006";
  this->VariableNames[1007] = "V1007";
  this->VariableNames[1008] = "V1008";
  this->VariableNames[1009] = "V1009";
  this->VariableNames[1010] = "V1010";
  this->VariableNames[1011] = "V1011";
  this->VariableNames[1012] = "V1012";
  this->VariableNames[1013] = "V1013";
  this->VariableNames[1014] = "V1014";
  this->VariableNames[1015] = "V1015";
  this->VariableNames[1016] = "V1016";
  this->VariableNames[1017] = "V1017";
  this->VariableNames[1018] = "V1018";
  this->VariableNames[1019] = "V1019";
  this->VariableNames[1020] = "V1020";
  this->VariableNames[1021] = "V1021";
  this->VariableNames[1022] = "V1022";
  this->VariableNames[1023] = "V1023";
  this->VariableNames[1024] = "V1024";
  this->VariableNames[1025] = "V1025";
  this->VariableNames[1026] = "V1026";
  this->VariableNames[1027] = "V1027";
  this->VariableNames[1028] = "V1028";
  this->VariableNames[1029] = "V1029";
  this->VariableNames[1030] = "V1030";
  this->VariableNames[1031] = "V1031";
  this->VariableNames[1032] = "V1032";
  this->VariableNames[1033] = "V1033";
  this->VariableNames[1034] = "V1034";
  this->VariableNames[1035] = "V1035";
  this->VariableNames[1036] = "V1036";
  this->VariableNames[1037] = "V1037";
  this->VariableNames[1038] = "V1038";
  this->VariableNames[1039] = "V1039";
  this->VariableNames[1040] = "V1040";
  this->VariableNames[1041] = "V1041";
  this->VariableNames[1042] = "V1042";
  this->VariableNames[1043] = "V1043";
  this->VariableNames[1044] = "V1044";
  this->VariableNames[1045] = "V1045";
  this->VariableNames[1046] = "V1046";
  this->VariableNames[1047] = "V1047";
  this->VariableNames[1048] = "V1048";
  this->VariableNames[1049] = "V1049";
  this->VariableNames[1050] = "V1050";
  this->VariableNames[1051] = "V1051";
  this->VariableNames[1052] = "V1052";
  this->VariableNames[1053] = "V1053";
  this->VariableNames[1054] = "V1054";
  this->VariableNames[1055] = "V1055";
  this->VariableNames[1056] = "V1056";
  this->VariableNames[1057] = "V1057";
  this->VariableNames[1058] = "V1058";
  this->VariableNames[1059] = "V1059";
  this->VariableNames[1060] = "V1060";
  this->VariableNames[1061] = "V1061";
  this->VariableNames[1062] = "V1062";
  this->VariableNames[1063] = "V1063";
  this->VariableNames[1064] = "V1064";
  this->VariableNames[1065] = "V1065";
  this->VariableNames[1066] = "V1066";
  this->VariableNames[1067] = "V1067";
  this->VariableNames[1068] = "V1068";
  this->VariableNames[1069] = "V1069";
  this->VariableNames[1070] = "V1070";
  this->VariableNames[1071] = "V1071";
  this->VariableNames[1072] = "V1072";
  this->VariableNames[1073] = "V1073";
  this->VariableNames[1074] = "V1074";
  this->VariableNames[1075] = "V1075";
  this->VariableNames[1076] = "V1076";
  this->VariableNames[1077] = "V1077";
  this->VariableNames[1078] = "V1078";
  this->VariableNames[1079] = "V1079";
  this->VariableNames[1080] = "V1080";
  this->VariableNames[1081] = "V1081";
  this->VariableNames[1082] = "V1082";
  this->VariableNames[1083] = "V1083";
  this->VariableNames[1084] = "V1084";
  this->VariableNames[1085] = "V1085";
  this->VariableNames[1086] = "V1086";
  this->VariableNames[1087] = "V1087";
  this->VariableNames[1088] = "V1088";
  this->VariableNames[1089] = "V1089";
  this->VariableNames[1090] = "V1090";
  this->VariableNames[1091] = "V1091";
  this->VariableNames[1092] = "V1092";
  this->VariableNames[1093] = "V1093";
  this->VariableNames[1094] = "V1094";
  this->VariableNames[1095] = "V1095";
  this->VariableNames[1096] = "V1096";
  this->VariableNames[1097] = "V1097";
  this->VariableNames[1098] = "V1098";
  this->VariableNames[1099] = "V1099";
  this->VariableNames[1100] = "V1100";
  this->VariableNames[1101] = "V1101";
  this->VariableNames[1102] = "V1102";
  this->VariableNames[1103] = "V1103";
  this->VariableNames[1104] = "V1104";
  this->VariableNames[1105] = "V1105";
  this->VariableNames[1106] = "V1106";
  this->VariableNames[1107] = "V1107";
  this->VariableNames[1108] = "V1108";
  this->VariableNames[1109] = "V1109";
  this->VariableNames[1110] = "V1110";
  this->VariableNames[1111] = "V1111";
  this->VariableNames[1112] = "V1112";
  this->VariableNames[1113] = "V1113";
  this->VariableNames[1114] = "V1114";
  this->VariableNames[1115] = "V1115";
  this->VariableNames[1116] = "V1116";
  this->VariableNames[1117] = "V1117";
  this->VariableNames[1118] = "V1118";
  this->VariableNames[1119] = "V1119";
  this->VariableNames[1120] = "V1120";
  this->VariableNames[1121] = "V1121";
  this->VariableNames[1122] = "V1122";
  this->VariableNames[1123] = "V1123";
  this->VariableNames[1124] = "V1124";
  this->VariableNames[1125] = "V1125";
  this->VariableNames[1126] = "V1126";
  this->VariableNames[1127] = "V1127";
  this->VariableNames[1128] = "V1128";
  this->VariableNames[1129] = "V1129";
  this->VariableNames[1130] = "V1130";
  this->VariableNames[1131] = "V1131";
  this->VariableNames[1132] = "V1132";
  this->VariableNames[1133] = "V1133";
  this->VariableNames[1134] = "V1134";
  this->VariableNames[1135] = "V1135";
  this->VariableNames[1136] = "V1136";
  this->VariableNames[1137] = "V1137";
  this->VariableNames[1138] = "V1138";
  this->VariableNames[1139] = "V1139";
  this->VariableNames[1140] = "V1140";
  this->VariableNames[1141] = "V1141";
  this->VariableNames[1142] = "V1142";
  this->VariableNames[1143] = "V1143";
  this->VariableNames[1144] = "V1144";
  this->VariableNames[1145] = "V1145";
  this->VariableNames[1146] = "V1146";
  this->VariableNames[1147] = "V1147";
  this->VariableNames[1148] = "V1148";
  this->VariableNames[1149] = "V1149";
  this->VariableNames[1150] = "V1150";
  this->VariableNames[1151] = "V1151";
  this->VariableNames[1152] = "V1152";
  this->VariableNames[1153] = "V1153";
  this->VariableNames[1154] = "V1154";
  this->VariableNames[1155] = "V1155";
  this->VariableNames[1156] = "V1156";
  this->VariableNames[1157] = "V1157";
  this->VariableNames[1158] = "V1158";
  this->VariableNames[1159] = "V1159";
  this->VariableNames[1160] = "V1160";
  this->VariableNames[1161] = "V1161";
  this->VariableNames[1162] = "V1162";
  this->VariableNames[1163] = "V1163";
  this->VariableNames[1164] = "V1164";
  this->VariableNames[1165] = "V1165";
  this->VariableNames[1166] = "V1166";
  this->VariableNames[1167] = "V1167";
  this->VariableNames[1168] = "V1168";
  this->VariableNames[1169] = "V1169";
  this->VariableNames[1170] = "V1170";
  this->VariableNames[1171] = "V1171";
  this->VariableNames[1172] = "V1172";
  this->VariableNames[1173] = "V1173";
  this->VariableNames[1174] = "V1174";
  this->VariableNames[1175] = "V1175";
  this->VariableNames[1176] = "V1176";
  this->VariableNames[1177] = "V1177";
  this->VariableNames[1178] = "V1178";
  this->VariableNames[1179] = "V1179";
  this->VariableNames[1180] = "V1180";
  this->VariableNames[1181] = "V1181";
  this->VariableNames[1182] = "V1182";
  this->VariableNames[1183] = "V1183";
  this->VariableNames[1184] = "V1184";
  this->VariableNames[1185] = "V1185";
  this->VariableNames[1186] = "V1186";
  this->VariableNames[1187] = "V1187";
  this->VariableNames[1188] = "V1188";
  this->VariableNames[1189] = "V1189";
  this->VariableNames[1190] = "V1190";
  this->VariableNames[1191] = "V1191";
  this->VariableNames[1192] = "V1192";
  this->VariableNames[1193] = "V1193";
  this->VariableNames[1194] = "V1194";
  this->VariableNames[1195] = "V1195";
  this->VariableNames[1196] = "V1196";
  this->VariableNames[1197] = "V1197";
  this->VariableNames[1198] = "V1198";
  this->VariableNames[1199] = "V1199";
  this->VariableNames[1200] = "V1200";
  this->VariableNames[1201] = "V1201";
  this->VariableNames[1202] = "V1202";
  this->VariableNames[1203] = "V1203";
  this->VariableNames[1204] = "V1204";
  this->VariableNames[1205] = "V1205";
  this->VariableNames[1206] = "V1206";
  this->VariableNames[1207] = "V1207";
  this->VariableNames[1208] = "V1208";
  this->VariableNames[1209] = "V1209";
  this->VariableNames[1210] = "V1210";
  this->VariableNames[1211] = "V1211";
  this->VariableNames[1212] = "V1212";
  this->VariableNames[1213] = "V1213";
  this->VariableNames[1214] = "V1214";
  this->VariableNames[1215] = "V1215";
  this->VariableNames[1216] = "V1216";
  this->VariableNames[1217] = "V1217";
  this->VariableNames[1218] = "V1218";
  this->VariableNames[1219] = "V1219";
  this->VariableNames[1220] = "V1220";
  this->VariableNames[1221] = "V1221";
  this->VariableNames[1222] = "V1222";
  this->VariableNames[1223] = "V1223";
  this->VariableNames[1224] = "V1224";
  this->VariableNames[1225] = "V1225";
  this->VariableNames[1226] = "V1226";
  this->VariableNames[1227] = "V1227";
  this->VariableNames[1228] = "V1228";
  this->VariableNames[1229] = "V1229";
  this->VariableNames[1230] = "V1230";
  this->VariableNames[1231] = "V1231";
  this->VariableNames[1232] = "V1232";
  this->VariableNames[1233] = "V1233";
  this->VariableNames[1234] = "V1234";
  this->VariableNames[1235] = "V1235";
  this->VariableNames[1236] = "V1236";
  this->VariableNames[1237] = "V1237";
  this->VariableNames[1238] = "V1238";
  this->VariableNames[1239] = "V1239";
  this->VariableNames[1240] = "V1240";
  this->VariableNames[1241] = "V1241";
  this->VariableNames[1242] = "V1242";
  this->VariableNames[1243] = "V1243";
  this->VariableNames[1244] = "V1244";
  this->VariableNames[1245] = "V1245";
  this->VariableNames[1246] = "V1246";
  this->VariableNames[1247] = "V1247";
  this->VariableNames[1248] = "V1248";
  this->VariableNames[1249] = "V1249";
  this->VariableNames[1250] = "V1250";
  this->VariableNames[1251] = "V1251";
  this->VariableNames[1252] = "V1252";
  this->VariableNames[1253] = "V1253";
  this->VariableNames[1254] = "V1254";
  this->VariableNames[1255] = "V1255";
  this->VariableNames[1256] = "V1256";
  this->VariableNames[1257] = "V1257";
  this->VariableNames[1258] = "V1258";
  this->VariableNames[1259] = "V1259";
  this->VariableNames[1260] = "V1260";
  this->VariableNames[1261] = "V1261";
  this->VariableNames[1262] = "V1262";
  this->VariableNames[1263] = "V1263";
  this->VariableNames[1264] = "V1264";
  this->VariableNames[1265] = "V1265";
  this->VariableNames[1266] = "V1266";
  this->VariableNames[1267] = "V1267";
  this->VariableNames[1268] = "V1268";
  this->VariableNames[1269] = "V1269";
  this->VariableNames[1270] = "V1270";
  this->VariableNames[1271] = "V1271";
  this->VariableNames[1272] = "V1272";
  this->VariableNames[1273] = "V1273";
  this->VariableNames[1274] = "V1274";
  this->VariableNames[1275] = "V1275";
  this->VariableNames[1276] = "V1276";
  this->VariableNames[1277] = "V1277";
  this->VariableNames[1278] = "V1278";
  this->VariableNames[1279] = "V1279";
  this->VariableNames[1280] = "V1280";
  this->VariableNames[1281] = "V1281";
  this->VariableNames[1282] = "V1282";
  this->VariableNames[1283] = "V1283";
  this->VariableNames[1284] = "V1284";
  this->VariableNames[1285] = "V1285";
  this->VariableNames[1286] = "V1286";
  this->VariableNames[1287] = "V1287";
  this->VariableNames[1288] = "V1288";
  this->VariableNames[1289] = "V1289";
  this->VariableNames[1290] = "V1290";
  this->VariableNames[1291] = "V1291";
  this->VariableNames[1292] = "V1292";
  this->VariableNames[1293] = "V1293";
  this->VariableNames[1294] = "V1294";
  this->VariableNames[1295] = "V1295";
  this->VariableNames[1296] = "V1296";
  this->VariableNames[1297] = "V1297";
  this->VariableNames[1298] = "V1298";
  this->VariableNames[1299] = "V1299";
  this->VariableNames[1300] = "V1300";
  this->VariableNames[1301] = "WSB";
  this->VariableNames[1302] = "WSN";
  this->VariableNames[1303] = "WSR";
  this->VariableNames[1304] = "WSB_M1";
  this->VariableNames[1305] = "WSB_M2";
  this->VariableNames[1306] = "WSN_M1";
  this->VariableNames[1307] = "WSN_M2";
  this->VariableNames[1308] = "WSR_M1";
  this->VariableNames[1309] = "WSR_M2";
  this->VariableNames[1310] = "MASGEN";
  this->VariableNames[1311] = "NUCRAT";
  this->VariableNames[1330] = "TEMPERATURE_M1";
  this->VariableNames[1331] = "TEMPERATURE_M2";
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCaseIndex(int ix)
{
  char b[5];

  if ( this->CaseFileBuffer[ix+2] == ' ')
    {
    b[0] = this->CaseFileBuffer[ix+1];
    b[1] = 0;
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( this->CaseFileBuffer[ix+3] == ' ')
    {
    b[0] = this->CaseFileBuffer[ix+1];
    b[1] = this->CaseFileBuffer[ix+2];
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( this->CaseFileBuffer[ix+4] == ' ')
    {
    b[0] = this->CaseFileBuffer[ix+1];
    b[1] = this->CaseFileBuffer[ix+2];
    b[2] = this->CaseFileBuffer[ix+3];
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if (this->CaseFileBuffer[ix+5] == ' ')
    {
    b[0] = this->CaseFileBuffer[ix+1];
    b[1] = this->CaseFileBuffer[ix+2];
    b[2] = this->CaseFileBuffer[ix+3];
    b[3] = this->CaseFileBuffer[ix+4];
    b[4] = 0;
    return atoi(b);
    }
  else
    {
    b[0] = this->CaseFileBuffer[ix+1];
    b[1] = this->CaseFileBuffer[ix+2];
    b[2] = this->CaseFileBuffer[ix+3];
    b[3] = this->CaseFileBuffer[ix+4];
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
      new_index = this->GetNothing(file_index);  // Section not used
      break;
    case 1:
      new_index = this->GetNothing(file_index);  // Section not used
      break;
    case 2:
      new_index = this->GetGridDimension(file_index);
      break;
    case 4:
      new_index = this->GetMachineConfiguration(file_index);
      break;
    case 10:
      new_index = this->GetNodesASCII(file_index);
      break;
    case 12:
      new_index = this->GetCellsASCII(file_index);
      break;
    case 13:
      new_index = this->GetFacesASCII(file_index);
      break;
    case 18:
      new_index = this->GetPeriodicShadowFacesASCII(file_index);
      break;
    case 33:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 37:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 38:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 39:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 40:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 41:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 45:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 54:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;
    case 58:
      new_index = this->GetCellTreeASCII(file_index);
      break;
    case 59:
      new_index = this->GetFaceTreeASCII(file_index);
      break;
    case 61:
      new_index = this->GetFaceParentsASCII(file_index);
      break;
    case 62:
      new_index = this->GetNCG1InformationASCII(file_index);
      break;
    case 63:
      new_index = this->GetNCG2InformationASCII(file_index);
      break;
    case 64:
      new_index = this->GetNoVariablesASCII(file_index);  // Section not used
      break;

    //
    // Single Precision
    //

    case 2010:
      new_index = this->GetNodesSinglePrecision(file_index);
      break;
    case 2012:
      new_index = this->GetCellsSinglePrecision(file_index);
      break;
    case 2013:
      new_index = this->GetFacesSinglePrecision(file_index);
      break;
    case 2018:
      new_index = this->GetPeriodicShadowFacesSinglePrecision(file_index);
      break;
    case 2033:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2033)");
      break;
    case 2037:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2037)");
      break;
    case 2038:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2038)");
      break;
    case 2039:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2039)");
      break;
    case 2040:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2040)");
      break;
    case 2041:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2041)");
      break;
    case 2045:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2045)");
      break;
    case 2058:
      new_index = this->GetCellTreeSinglePrecision(file_index);
      break;
    case 2059:
      new_index = this->GetFaceTreeSinglePrecision(file_index);
      break;
    case 2061:
      new_index = this->GetFaceParentsSinglePrecision(file_index);
      break;
    case 2062:
      new_index = this->GetNCG1InformationSinglePrecision(file_index);
      break;
    case 2063:
      new_index = this->GetNCG2InformationSinglePrecision(file_index);
      break;
    case 2064:
      new_index = this->GetNoVariablesSinglePrecision(file_index, "2064)");
      break;

    case 3010:
      new_index = this->GetNodesDoublePrecision(file_index);
      break;
    case 3012:
      new_index = this->GetCellsDoublePrecision(file_index);
      break;
    case 3013:
      new_index = this->GetFacesDoublePrecision(file_index);
      break;
    case 3018:
      new_index = this->GetPeriodicShadowFacesDoublePrecision(file_index);
      break;
    case 3033:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3033)");
      break;
    case 3037:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3037)");
      break;
    case 3038:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3038)");
      break;
    case 3039:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3039)");
      break;
    case 3040:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3040)");
      break;
    case 3041:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3041)");
      break;
    case 3045:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3045)");
      break;
    case 3058:
      new_index = this->GetCellTreeDoublePrecision(file_index);
      break;
    case 3059:
      new_index = this->GetFaceTreeDoublePrecision(file_index);
      break;
    case 3061:
      new_index = this->GetFaceParentsDoublePrecision(file_index);
      break;
    case 3062:
      new_index = this->GetNCG1InformationDoublePrecision(file_index);
      break;
    case 3063:
      new_index = this->GetNCG2InformationDoublePrecision(file_index);
      break;
    case 3064:
      new_index = this->GetNoVariablesDoublePrecision(file_index, "3064)");
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

  if ( this->DataFileBuffer[ix+2] == ' ')
    {
    b[0] = this->DataFileBuffer[ix+1];
    b[1] = 0;
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( this->DataFileBuffer[ix+3] == ' ')
    {
    b[0] = this->DataFileBuffer[ix+1];
    b[1] = this->DataFileBuffer[ix+2];
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( this->DataFileBuffer[ix+4] == ' ')
    {
    b[0] = this->DataFileBuffer[ix+1];
    b[1] = this->DataFileBuffer[ix+2];
    b[2] = this->DataFileBuffer[ix+3];
    b[3] = 0;
    b[4] = 0;
    return atoi(b);
    }
  else if ( this->DataFileBuffer[ix+5] == ' ')
    {
    b[0] = this->DataFileBuffer[ix+1];
    b[1] = this->DataFileBuffer[ix+2];
    b[2] = this->DataFileBuffer[ix+3];
    b[3] = this->DataFileBuffer[ix+4];
    b[4] = 0;
    return atoi(b);
    }
  else
    {
    b[0] = this->DataFileBuffer[ix+1];
    b[1] = this->DataFileBuffer[ix+2];
    b[2] = this->DataFileBuffer[ix+3];
    b[3] = this->DataFileBuffer[ix+4];
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
      new_index = this->GetDataNothing(file_index);
      break;
    case 1:
      new_index = this->GetDataNothing(file_index);
      break;
    case 2:
      new_index = this->GetDataGridDimension(file_index);
      break;
    case 4:
      new_index = this->GetNoData(file_index);
      break;
    case 33:
      new_index = this->GetNoData(file_index);
      break;
    case 37:
      new_index = this->GetNoData(file_index);
      break;
    case 38:
      new_index = this->GetNoData(file_index);
      break;
    case 50:
      new_index = this->GetDataUnknownASCII(file_index);
      break;
    case 64:
      new_index = this->GetNoData(file_index);
      break;
    case 300:
      new_index = this->GetDataASCII(file_index);
      break;
    case 301:
      new_index = this->GetNoData(file_index);
      break;
    case 302:
      new_index = this->GetNoData(file_index);
      break;
    case 303:
      new_index = this->GetNoData(file_index);
      break;
    case 313:
      new_index = this->GetNoData(file_index);
      break;
    case 2300:
      new_index = this->GetDataSinglePrecision(file_index);
      break;
    case 2301:
      new_index = this->SkipUnknownSinglePrecisionData(file_index, "2301)");
      break;
    case 2302:
      new_index = this->SkipUnknownSinglePrecisionData(file_index, "2302)");
      break;
    case 2303:
      new_index = this->SkipUnknownSinglePrecisionData(file_index, "2303)");
      break;
    case 2313:
      new_index = this->SkipUnknownSinglePrecisionData(file_index, "2313)");
      break;
    case 3300:
      new_index = this->GetDataDoublePrecision(file_index);
      break;
    case 3301:
      new_index = this->SkipUnknownDoublePrecisionData(file_index, "3301)");
      break;
    case 3302:
      new_index = this->SkipUnknownDoublePrecisionData(file_index, "3302)");
      break;
    case 3313:
      new_index = this->SkipUnknownDoublePrecisionData(file_index, "3313)");
      break;
    default:
      cout << " Unknown Index " << task << endl;
      exit(1);
      break;
    }
  return new_index;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNothing(int ix)
{
  return this->GoToNextRightParen(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetMachineConfiguration(int ix)
{
  char buf[120];
  int j = ix+1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );
  j = this->GoToNextRightParen(j)+1;

  int a, b, c, d, e, f, g, h, m, n, o;
  sscanf( buf, " %d %d %d %d %d %d %d %d %d %d %d", 
    &a, &b, &c, &d, &e, &f, &g, &h, &m, &n, &o );

  if ( a == 60 )
    {
    this->LittleEndianFlag = 1;
    }
  else
    {
    this->LittleEndianFlag = 0;
    }

  return this->GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNoVariablesASCII(int ix)
{
  return this->GoToNextSectionASCII(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsASCII(int ix)
{
  char buf[120];
  int j = ix+1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );
  j = this->GoToNextRightParen(j)+1;

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if ( zi != 0)
    {
    this->CellZones->InsertValue(NumberOfCellZones, zi);
    this->NumberOfCellZones++;
    }

  if ( zi == 0) 
    {
    this->NumberOfCells = li;
    this->CellParentFlags.resize( NumberOfCells+1, false ); // mccdo
    }
  else
    {
    if ( et == 0)
      {
      this->GetMixedCellTypes( j, fi, li); 
      }
    else
      {
      for (int i = fi; i <= li; i++)
        {
        this->CellTypes->InsertValue(i, et);
        }
      }
    }
  return this->GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if (zi == 0)
    {
    this->NumberOfFaces = li;
    this->InterfaceFaceChildFlags.resize( NumberOfFaces+1, false ); //mccdo
    this->FaceParentFlags.resize( NumberOfFaces+1, false ); //mccdo
    }
  else
    {
    j = this->GoToNextLeftParen(j)+1;
    j = this->GoToNextEOL(j) +1;
    int n0, n1, n2, n3;
    int c0, c1;
    int type;
    for (int k = fi; k <= li; k++)
      {
      this->GetStringToNextRightParenOrEOL( j, buf );
      if ( et == 0 )
        {
        if ( buf[0] == 2)
          {
          sscanf( buf, " %x %x %x %x %x ", &type, &n0 , &n1, &c0, &c1 );
          this->FaceTypes->InsertValue(k,type);
          this->FaceNodes->InsertComponent(k,0,n0-1);
          this->FaceNodes->InsertComponent(k,1,n1-1);
          this->FaceNodes->InsertComponent(k,2,0);
          this->FaceNodes->InsertComponent(k,3,0);
          this->FaceCells->InsertComponent(k,0,c0);
          this->FaceCells->InsertComponent(k,1,c1);
          }
        else if ( buf[1] == 3)
          {
          sscanf( buf, " %x %x %x %x %x %x ", &type , &n0, &n1, 
            &n2, &c0, &c1 );
          this->FaceTypes->InsertValue(k,type);
          this->FaceNodes->InsertComponent(k,0,n0-1);
          this->FaceNodes->InsertComponent(k,1,n1-1);
          this->FaceNodes->InsertComponent(k,2,n2-1);
          this->FaceNodes->InsertComponent(k,3,0);
          this->FaceCells->InsertComponent(k,0,c0);
          this->FaceCells->InsertComponent(k,1,c1);
          }
        else
          {
          sscanf( buf, " %x %x %x %x %x %x %x ", &type, &n0 , &n1,
            &n2, &n3, &c0, &c1 );
          this->FaceTypes->InsertValue(k,type);
          this->FaceNodes->InsertComponent(k,0,n0-1);
          this->FaceNodes->InsertComponent(k,1,n1-1);
          this->FaceNodes->InsertComponent(k,2,n2-1);
          this->FaceNodes->InsertComponent(k,3,n3-1);
          this->FaceCells->InsertComponent(k,0,c0);
          this->FaceCells->InsertComponent(k,1,c1);
          }
        }
      else if ( et == 2)
        {
        sscanf( buf, " %x %x %x %x ", &n0 , &n1, &c0, &c1 );
        this->FaceTypes->InsertValue(k,2);
        this->FaceNodes->InsertComponent(k,0,n0-1);
        this->FaceNodes->InsertComponent(k,1,n1-1);
        this->FaceNodes->InsertComponent(k,2,0);
        this->FaceNodes->InsertComponent(k,3,0);
        this->FaceCells->InsertComponent(k,0,c0);
        this->FaceCells->InsertComponent(k,1,c1);
        }
      else if ( et == 3)
        {
        sscanf( buf, " %x %x %x %x %x ", &n0 , &n1, &n2, &c0, &c1 );
        this->FaceTypes->InsertValue(k,3);
        this->FaceNodes->InsertComponent(k,0,n0-1);
        this->FaceNodes->InsertComponent(k,1,n1-1);
        this->FaceNodes->InsertComponent(k,2,n2-1);
        this->FaceNodes->InsertComponent(k,3,0);
        this->FaceCells->InsertComponent(k,0,c0);
        this->FaceCells->InsertComponent(k,1,c1);
        }
      else
        {
        sscanf( buf, " %x %x %x %x %x %x ", &n0 , &n1, &n2, &n3, &c0, &c1 );
        this->FaceTypes->InsertValue(k,4);
        this->FaceNodes->InsertComponent(k,0,n0-1);
        this->FaceNodes->InsertComponent(k,1,n1-1);
        this->FaceNodes->InsertComponent(k,2,n2-1);
        this->FaceNodes->InsertComponent(k,3,n3-1);
        this->FaceCells->InsertComponent(k,0,c0);
        this->FaceCells->InsertComponent(k,1,c1);
        }
      j = this->GoToNextEOL(j) +1;
      }
    }
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  this->Points->InsertPoint(0, 0.0 , 0.0 , 0.0);

  if (zi == 0)
    {
    this->NumberOfNodes = li;
    }
  else
    {
    j = this->GoToNextLeftParen(j)+1;
    j = this->GoToNextEOL(j) +1;
    float x,y,z;
    for (int k = fi; k <= li; k++)
      {
      this->GetStringToNextRightParenOrEOL( j, buf );
      if ( nd == 2)
        {
        sscanf( buf, " %f %f ", &x , &y );
        this->Points->InsertPoint(k-1, x, y, 0.0);
        }
      else
        {
        sscanf( buf, " %f %f %f", &x , &y, &z );
        this->Points->InsertPoint(k-1, x, y, z);
        }

      j = this->GoToNextEOL(j) +1;
      }
    }
  j = this->GoToNextRightParen(j)+1;
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsASCII(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);

  j = this->GoToNextLeftParen(j)+1;
  j = this->GoToNextASCIIHexDigit(j);

  for (int k=face_id0;k<=face_id1;k++)
    {
    this->GetStringToNextRightParenOrEOL( j, buf );
    int pid0, pid1;
    sscanf( buf, " %x %x ", &pid0 , &pid1 );

    // mccdo this->FaceParents->InsertComponent(k, 0, pid0);
    // mccdo this->FaceParents->InsertComponent(k, 1, pid1);
    // mccdo this->FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    InterfaceFaceChildFlags[ k ] = true; // mccdo
    this->NumberOfFaceParentChildren++;

    j = this->GoToNextEOL(j) +1;
    }

  if (face_id1 >= this->NumberOfFaceParents)
    {
    this->NumberOfFaceParents = face_id1;
    }

  return this->GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationASCII(int ix)
{
  // Face Information
  char buf[120];
  int j = ix + 1;

  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );

  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);

  j = this->GoToNextLeftParen(j)+1;
  j = this->GoToNextASCIIHexDigit(j);

  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    this->GetStringToNextRightParenOrEOL( j, buf );
    int child, parent;
    sscanf( buf, " %x %x ", &child , &parent );
    // mccdo this->NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    // mccdo this->NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    this->NCGFaceChildFlags.insert( child ); // mccdo
    j = this->GoToNextEOL(j) +1;
    this->NumberOfNCGFaces++;
    }

  this->NumberOfNCGFaceHeaders++;
  return this->GoToNextRightParen(j)+1;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationASCII(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);

  j = this->GoToNextLeftParen(j)+1;
  j = this->GoToNextASCIIHexDigit(j);

  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    this->GetStringToNextRightParenOrEOL( j, buf );
    float x,y,z;
    int NodeId;
    if (this->GridDimension == 3)
      {
      sscanf( buf, " %d %f %f %f ", &NodeId, &x , &y, &z );
      }
    else
      {
      sscanf( buf, " %d %f %f ", &NodeId, &x , &y );
      z = 0;
      }


    j = this->GoToNextEOL(j) +1;
    this->NumberOfNCGNodes++;
    }

  this->NumberOfNCGNodeHeaders++;
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesASCII(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = this->GoToNextLeftParen(j)+1;
  j = this->GoToNextASCIIHexDigit(j);

  int psf0, psf1;
  for (int k = fi; k <= li; k++)
    {
    this->GetStringToNextRightParenOrEOL( j, buf );
    sscanf( buf, " %x %x ", &psf0 , &psf1 );
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    j = this->GoToNextEOL(j) +1;
    }

  if ( li >= this->NumberOfPeriodicShadowFaces)
    {
    this->NumberOfPeriodicShadowFaces = li;
    }

  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeASCII(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  // mccdo this->CellTreeParentCellId0->InsertValue(this->NumberOfCellTrees, fid0);
  // mccdo this->CellTreeParentCellId1->InsertValue(this->NumberOfCellTrees, fid1);
  // mccdo
  for (int k = fid0; k <= fid1; k++)
    {
    this->CellParentFlags.at( k ) = true;
    }
  // mccdo
  j = this->GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetAsciiInteger(j);
    j = this->GoPastAsciiInteger(j);
    // mccdo this->CellTreesNumberOfKids->InsertValue(this->NumberOfCellTreeParents, 
    // mccdo   NumberOfKids);
    // mccdo this->CellTreesKidsIndex->InsertValue(this->NumberOfCellTreeParents, 
    // mccdo   this->NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      //int Kid = this->GetAsciiInteger(j);
      j = this->GoPastAsciiInteger(j);
      // mccdo this->CellTreesKids->InsertValue(this->NumberOfCellTreeKids, Kid);
      this->NumberOfCellTreeKids++;
      }
    this->NumberOfCellTreeParents++;
    }

  this->NumberOfCellTrees++;
  return this->GoToNextSectionASCII(j);
}

int vtkFLUENTReader::GetFaceTreeASCII(int ix)
{
  char buf[120];
  int j = ix + 1;

  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  // mccdo this->FaceTreeParentFaceId0->InsertValue(this->NumberOfFaceTrees, fid0);
  // mccdo this->FaceTreeParentFaceId1->InsertValue(this->NumberOfFaceTrees, fid1);
  // mccdo
  static int index = 0;
  for(int k = fid0; k <= fid1; k++)
    {
    this->FaceTreeParentTable->InsertValue(k, index);
    index++;
    }

  int startFace = fid0;
  int endFace = fid1;
  for(int k = startFace; k <= endFace; k++)
    {
    this->FaceParentFlags.at( k ) = true;
    }
  //mccdo
  j = this->GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetAsciiInteger(j);
    j = this->GoPastAsciiInteger(j);
    this->FaceTreesNumberOfKids->InsertValue(this->NumberOfFaceTreeParents, 
      NumberOfKids);
    this->FaceTreesKidsIndex->InsertValue(this->NumberOfFaceTreeParents, 
      this->NumberOfFaceTreeKids);
    for(int i = 0; i < NumberOfKids; i++)
      {
      int Kid = this->GetAsciiInteger(j);
      j = this->GoPastAsciiInteger(j);
      this->FaceTreesKids->InsertValue(this->NumberOfFaceTreeKids, Kid);
      this->NumberOfFaceTreeKids++;
      }

    this->NumberOfFaceTreeParents++;
    }

  this->NumberOfFaceTrees++;
  return this->GoToNextSectionASCII(j);
}

//-----------------------------------------------------------------------------

int vtkFLUENTReader::GetNoVariablesSinglePrecision(int ix, char buf[])
{
  return this->GoToNextSectionSinglePrecision( ix, buf);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );

  if ( zi != 0)
    {
    this->CellZones->InsertValue(this->NumberOfCellZones, zi);
    this->NumberOfCellZones++;
    }

  if ( et != 0)
    {
    for (int i = fi; i <= li; i++)
      {
      this->CellTypes->InsertValue(i, et);
      }
    }
  else
    { // Mixed Cells
    j = this->GoToNextLeftParen(j)+1;
    for (int i = fi; i <= li; i++)
      {
      this->CellTypes->InsertValue(i, this->GetBinaryInteger(j));
      j = j + 4;
      }
    }

  j++;
  return this->GoToNextSectionSinglePrecision( j, "2012)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  j = this->GoToNextLeftParen(j)+1;

  if ( et == 2)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2, 0);
      this->FaceNodes->InsertComponent(i,3, 0);
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if ( et == 3)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,3, 0);
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if ( et == 4)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,3,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else
    { // Mixed Faces
    for (int i = fi; i <= li; i++)
      {
      int ft = this->GetBinaryInteger(j);
      j = j + 4;
      this->FaceTypes->InsertValue(i, ft);
      if ( ft == 2)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2, 0);
        this->FaceNodes->InsertComponent(i,3, 0);
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 3)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,3, 0);
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 4)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,3,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      }
    }
  return this->GoToNextSectionSinglePrecision( j, "2013)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  this->Points->InsertPoint(0, 0.0 , 0.0 , 0.0);
  j = this->GoToNextLeftParen(j)+1;

  float x,y,z;
  for(int k = fi; k <= li; k++)
    {
    if ( nd == 2)
      {
      x = this->GetBinaryFloat(j);
      j = j+4;
      y = this->GetBinaryFloat(j);
      j = j+4;
      this->Points->InsertPoint(k-1, x, y, 0);
      }
    else
      {
      x = this->GetBinaryFloat(j);
      j = j+4;
      y = this->GetBinaryFloat(j);
      j = j+4;
      z = this->GetBinaryFloat(j);
      j = j+4;
      this->Points->InsertPoint(k-1, x, y, z);
      }
    }

  return this->GoToNextSectionSinglePrecision( j, "2010)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);
  j = this->GoToNextLeftParen(j)+1;

  int pid0, pid1;
  for (int k = face_id0; k <= face_id1; k++)
    {
    pid0 = this->GetBinaryInteger(j);
    j = j + 4;
    pid1 = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->FaceParents->InsertComponent(k, 0, pid0);
    // mccdo this->FaceParents->InsertComponent(k, 1, pid1);
    // mccdo this->FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    InterfaceFaceChildFlags[ k ] = true; // mccdo
    this->NumberOfFaceParentChildren++;
    }

  if ( face_id1 >= this->NumberOfFaceParents)
    {
    this->NumberOfFaceParents = face_id1;
    }

  return this->GoToNextSectionSinglePrecision( j, "2061)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);

  j = this->GoToNextLeftParen(j)+1;
  int child,parent;
  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    child = this->GetBinaryInteger(j);
    j = j + 4;
    parent = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    // mccdo this->NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    this->NCGFaceChildFlags.insert( child ); // mccdo
    this->NumberOfNCGFaces++;
    }

  this->NumberOfNCGFaceHeaders++;
  return this->GoToNextSectionSinglePrecision( j, "2062)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationSinglePrecision(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);

  j = this->GoToNextLeftParen(j)+1;

  float x,y,z;
  int NodeId;
  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    NodeId = this->GetBinaryInteger(j);
    j = j + 4;
    x = this->GetBinaryFloat(j);
    j = j + 4;
    y = this->GetBinaryFloat(j);
    j = j + 4;
    if ( this->GridDimension == 3)
      {
      z = this->GetBinaryFloat(j);
      j = j + 4;
      }
    else
      {
      z = 0.0;
      }

    this->NumberOfNCGNodes++;
    }

  this->NumberOfNCGNodeHeaders++;
  return this->GoToNextSectionSinglePrecision( j, "2063)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = this->GoToNextLeftParen(j)+1;

  int psf0, psf1;
  for (int k = fi; k <= li; k++)
    {
    psf0 = this->GetBinaryInteger(j);
    j = j + 4;
    psf1 = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    }

  if ( li >= this->NumberOfPeriodicShadowFaces)
    {
    this->NumberOfPeriodicShadowFaces = li;
    }

  return this->GoToNextSectionSinglePrecision( j, "2018)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  // mccdo this->CellTreeParentCellId0->InsertValue(this->NumberOfCellTrees, fid0);
  // mccdo this->CellTreeParentCellId1->InsertValue(this->NumberOfCellTrees, fid1);
  // mccdo
  for (int k = fid0; k <= fid1; k++)
    {
    this->CellParentFlags.at( k ) = true;
    }
  // mccdo
  j = this->GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->CellTreesNumberOfKids->InsertValue(this->NumberOfCellTreeParents,
    // mccdo   NumberOfKids);
    // mccdo this->CellTreesKidsIndex->InsertValue(this->NumberOfCellTreeParents, 
    // mccdo   NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      //int Kid = this->GetBinaryInteger(j);
      j = j + 4;
      // mccdo this->CellTreesKids->InsertValue(this->NumberOfCellTreeKids, Kid);
      this->NumberOfCellTreeKids++;
      }

    this->NumberOfCellTreeParents++;
    }

  this->NumberOfCellTrees++;
  return this->GoToNextSectionSinglePrecision( j, "2058)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceTreeSinglePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);

  // mccdo this->FaceTreeParentFaceId0->InsertValue(this->NumberOfFaceTrees, fid0);
  // mccdo this->FaceTreeParentFaceId1->InsertValue(this->NumberOfFaceTrees, fid1);
  // mccdo
  static int index = 0;
  for(int k = fid0; k <= fid1; k++)
    {
    this->FaceTreeParentTable->InsertValue(k, index);
    index++;
    }

  int startFace = fid0;
  int endFace = fid1;
  for(int k = startFace; k <= endFace; k++)
    {
    this->FaceParentFlags.at( k ) = true;
    }
  //mccdo
  j = this->GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetBinaryInteger(j);
    j = j + 4;
    this->FaceTreesNumberOfKids->InsertValue(this->NumberOfFaceTreeParents, 
      NumberOfKids);
    this->FaceTreesKidsIndex->InsertValue(this->NumberOfFaceTreeParents, 
      NumberOfFaceTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = this->GetBinaryInteger(j);
      j = j + 4;
      this->FaceTreesKids->InsertValue(this->NumberOfFaceTreeKids, Kid);
      this->NumberOfFaceTreeKids++;
      }
    this->NumberOfFaceTreeParents++;
    }
  this->NumberOfFaceTrees++;
  return this->GoToNextSectionSinglePrecision( j, "2059)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNoVariablesDoublePrecision(int ix, char buf[])
{
  return this->GoToNextSectionDoublePrecision( ix, buf);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellsDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  if ( zi != 0)
    {
    this->CellZones->InsertValue(this->NumberOfCellZones, zi);
    this->NumberOfCellZones++;
    }

  if ( et != 0)
    {
    for ( int i = fi; i <= li; i++)
      {
      this->CellTypes->InsertValue(i, et);
      }
    }
  else
    { // Mixed Cells
    j = this->GoToNextLeftParen(j)+1;
    for (int i = fi; i <= li; i++)
      {
      this->CellTypes->InsertValue(i, this->GetBinaryInteger(j));
      j = j + 4;
      }
    }

  j++;
  return this->GoToNextSectionDoublePrecision( j, "3012)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFacesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, et;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &et );
  j = this->GoToNextLeftParen(j)+1;

  if (et == 2)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2, 0);
      this->FaceNodes->InsertComponent(i,3, 0);
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if (et == 3)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,3, 0);
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else if (et == 4)
    {
    for (int i = fi; i <= li; i++)
      {
      this->FaceTypes->InsertValue(i, et);
      this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceNodes->InsertComponent(i,3,this->GetBinaryInteger(j)-1);
      j = j + 4;
      this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
      j = j + 4;
      this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
      j = j + 4;
      }
    }
  else
    { // Mixed Faces
    for (int i = fi; i <= li; i++)
      {
      int ft = this->GetBinaryInteger(j);
      j = j + 4;
      this->FaceTypes->InsertValue(i, ft);
      if ( ft == 2)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2, 0);
        this->FaceNodes->InsertComponent(i,3, 0);
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 3)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,3, 0);
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      else if ( ft == 4)
        {
        this->FaceNodes->InsertComponent(i,0,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,1,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,2,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceNodes->InsertComponent(i,3,this->GetBinaryInteger(j)-1);
        j = j + 4;
        this->FaceCells->InsertComponent(i,0,this->GetBinaryInteger(j));
        j = j + 4;
        this->FaceCells->InsertComponent(i,1,this->GetBinaryInteger(j));
        j = j + 4;
        }
      }
    }
  return this->GoToNextSectionDoublePrecision( j, "3013)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNodesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int zi, fi, li, ty, nd;
  sscanf( buf, " %x %x %x %x %x", &zi, &fi, &li, &ty, &nd );

  this->Points->InsertPoint(0, 0.0 , 0.0 , 0.0);
  j = this->GoToNextLeftParen(j)+1;

  float x,y,z;
  for (int k = fi; k <= li; k++)
    {
    if ( nd == 2)
      {
      x = this->GetBinaryDouble(j);
      j = j+8;
      y = this->GetBinaryDouble(j);
      j = j+8;
      this->Points->InsertPoint(k-1, x, y, 0);
      }
    else
      {
      x = this->GetBinaryDouble(j);
      j = j+8;
      y = this->GetBinaryDouble(j);
      j = j+8;
      z = this->GetBinaryDouble(j);
      j = j+8;
      this->Points->InsertPoint(k-1, x, y, z);
      }
    }
  return this->GoToNextSectionSinglePrecision( j, "3010)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceParentsDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int face_id0, face_id1;
  sscanf( buf, " %x %x", &face_id0, &face_id1);
  j = this->GoToNextLeftParen(j)+1;

  int pid0, pid1;
  for (int k = face_id0; k <= face_id1; k++)
    {
    pid0 = this->GetBinaryInteger(j);
    j = j + 4;
    pid1 = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->FaceParents->InsertComponent(k, 0, pid0);
    // mccdo this->FaceParents->InsertComponent(k, 1, pid1);
    // mccdo this->FaceParentsChildren->InsertValue(NumberOfFaceParentChildren, k);
    InterfaceFaceChildFlags[ k ] = true; // mccdo
    this->NumberOfFaceParentChildren++;
    }

  if (face_id1 >= this->NumberOfFaceParents)
    {
    this->NumberOfFaceParents = face_id1;
    }

  return this->GoToNextSectionDoublePrecision( j, "3061)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG1InformationDoublePrecision(int ix)
{
  // Face Information
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;

  this->GetStringToNextRightParen( j, buf );
  int KidId, ParentId, NumberOfFacesNCG;
  sscanf( buf, " %d %d %d", &KidId, &ParentId, &NumberOfFacesNCG);

  j = this->GoToNextLeftParen(j)+1;
  int child,parent;
  for (int k = 0; k < NumberOfFacesNCG; k++)
    {
    child = this->GetBinaryInteger(j);
    j = j + 4;
    parent = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->NCGFaceChild->InsertValue(NumberOfNCGFaces, child);
    // mccdo this->NCGFaceParent->InsertValue(NumberOfNCGFaces, parent);
    this->NCGFaceChildFlags.insert( child ); //mccdo 
    this->NumberOfNCGFaces++;
    }

  this->NumberOfNCGFaceHeaders++;
  return this->GoToNextSectionDoublePrecision( j, "3062)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNCG2InformationDoublePrecision(int ix)
{
  // Node Information
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int ZoneId, NumberOfNodesNCG;
  sscanf( buf, " %d %d", &ZoneId, &NumberOfNodesNCG);

  j = this->GoToNextLeftParen(j)+1;

  float x,y,z;
  int NodeId;
  for (int k = 0; k < NumberOfNodesNCG; k++)
    {
    NodeId = this->GetBinaryInteger(j);
    j = j + 4;
    x = this->GetBinaryDouble(j);
    j = j + 8;
    y = this->GetBinaryDouble(j);
    j = j + 8;
    if (this->GridDimension == 3)
      {
      z = this->GetBinaryDouble(j);
      j = j + 8;
      }
    else
      {
      z = 0.0;
      }

    this->NumberOfNCGNodes++;
    }

  this->NumberOfNCGNodeHeaders++;
  return this->GoToNextSectionDoublePrecision( j, "3063)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetPeriodicShadowFacesDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );
  int fi, li, pz, sz;
  sscanf( buf, " %x %x %x %x", &fi, &li, &pz, &sz);
  j = this->GoToNextLeftParen(j)+1;

  int psf0, psf1;
  for (int k = fi;k <= li; k++)
    {
    psf0 = this->GetBinaryInteger(j);
    j = j + 4;
    psf1 = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 0, psf0);
    // mccdo this->PeriodicShadowFaces->InsertComponent(k, 1, psf1);
    }

  if ( li >= this->NumberOfPeriodicShadowFaces)
    {
    this->NumberOfPeriodicShadowFaces = li;
    }

  return this->GoToNextSectionDoublePrecision( j, "3018)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetCellTreeDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);
  // mccdo this->CellTreeParentCellId0->InsertValue(this->NumberOfCellTrees, fid0);
  // mccdo this->CellTreeParentCellId1->InsertValue(this->NumberOfCellTrees, fid1);
  // mccdo
  for (int k = fid0; k <= fid1; k++)
    {
    this->CellParentFlags.at( k ) = true;
    }
  // mccdo

  j = this->GoToNextLeftParen(j)+1;
  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetBinaryInteger(j);
    j = j + 4;
    // mccdo this->CellTreesNumberOfKids->InsertValue(this->NumberOfCellTreeParents,
    // mccdo   NumberOfKids);
    // mccdo this->CellTreesKidsIndex->InsertValue(this->NumberOfCellTreeParents, 
    // mccdo   NumberOfCellTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      //int Kid = this->GetBinaryInteger(j);
      j = j + 4;
      // mccdo this->CellTreesKids->InsertValue(this->NumberOfCellTreeKids, Kid);
      this->NumberOfCellTreeKids++;
      }
    this->NumberOfCellTreeParents++;
    }
  this->NumberOfCellTrees++;
  return this->GoToNextSectionDoublePrecision( j, "3058)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetFaceTreeDoublePrecision(int ix)
{
  char buf[120];
  int j = ix + 1;
  j = this->GoToNextLeftParen(j)+1;
  this->GetStringToNextRightParen( j, buf );

  int fid0, fid1, pzid, czid;
  sscanf( buf, " %x %x %x %x", &fid0, &fid1, &pzid, &czid);
  // mccdo this->FaceTreeParentFaceId0->InsertValue(this->NumberOfFaceTrees, fid0);
  // mccdo this->FaceTreeParentFaceId1->InsertValue(this->NumberOfFaceTrees, fid1);
  // mccdo
  static int index = 0;
  for(int k = fid0; k <= fid1; k++)
    {
    this->FaceTreeParentTable->InsertValue(k, index);
    index++;
    }

  int startFace = fid0;
  int endFace = fid1;
  for(int k = startFace; k <= endFace; k++)
    {
    this->FaceParentFlags.at( k ) = true;
    }
  //mccdo
  j = GoToNextLeftParen(j)+1;

  for (int k = fid0; k <= fid1; k++)
    {
    int NumberOfKids = this->GetBinaryInteger(j);
    j = j + 4;
    this->FaceTreesNumberOfKids->InsertValue(this->NumberOfFaceTreeParents,
      NumberOfKids);
    this->FaceTreesKidsIndex->InsertValue(this->NumberOfFaceTreeParents, 
      this->NumberOfFaceTreeKids);
    for (int i = 0; i < NumberOfKids; i++)
      {
      int Kid = this->GetBinaryInteger(j);
      j = j + 4;
      this->FaceTreesKids->InsertValue(this->NumberOfFaceTreeKids, Kid);
      this->NumberOfFaceTreeKids++;
      }
    this->NumberOfFaceTreeParents++;
    }
  this->NumberOfFaceTrees++;
  return this->GoToNextSectionDoublePrecision( j, "3059)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetGridDimension(int ix)
{
  char b2[2];
  b2[0] = this->CaseFileBuffer[ix+3];
  b2[1] = 0;
  this->GridDimension = atoi(b2);

  return this->GoToNextRightParen(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataNothing(int ix)
{
  return this->GoToNextRightParenData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetNoData(int ix)
{
  return this->GoToNextSectionASCIIData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataGridDimension(int ix)
{
  char b2[2];
  b2[0] = this->DataFileBuffer[ix+3];
  b2[1] = 0;
  this->GridDimension = atoi(b2);

  return this->GoToNextRightParenData(ix);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextRightParenData(int ix)
{
  int i = ix;
  while ( this->DataFileBuffer[i] != ')' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextLeftParenData(int ix)
{
  int i = ix;
  while ( this->DataFileBuffer[i] != '(' )
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
    if ( this->DataFileBuffer[i] == '(')
      {
      level++;
      }
    if (this->DataFileBuffer[i] == ')')
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
  while ( !((this->DataFileBuffer[i] == buf[0]) 
    && (this->DataFileBuffer[i+1] == buf[1]) 
    && (this->DataFileBuffer[i+2] == buf[2]) 
    && (this->DataFileBuffer[i+3] == buf[3]) 
    && (this->DataFileBuffer[i+4] == buf[4])))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionDoublePrecisionData(int ix, char buf[])
{
  int i = ix + 1;
  while( !((this->DataFileBuffer[i] == buf[0]) 
    && (this->DataFileBuffer[i+1] == buf[1]) 
    && (this->DataFileBuffer[i+2] == buf[2]) 
    && (this->DataFileBuffer[i+3] == buf[3]) 
    && (this->DataFileBuffer[i+4] == buf[4])))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataASCII(int ix)
{
  char buf[120];
  int index = -1; // mccdo
  int j = ix + 1;
  float x;
  j = this->GoToNextLeftParenData(j)+1;
  this->GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = this->GoToNextLeftParenData(j)+1;

  if (this->DataPass == 1)
    {
    if (this->IsCellZoneId(zid))
      {
      if (this->IsNewVariable(ssid))
        {
        this->VariableIds->InsertValue(NumberOfVariables, ssid);
        this->VariableSizes->InsertValue(NumberOfVariables, size);
        this->NumberOfVariables++;
        }
      }
    }
  else 
    {
    if (this->IsCellZoneId(zid))
      {
      index = this->GetVariableIndex(ssid); // mccdo
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          this->GetStringToNextRightParenOrEOLData( j, buf );
          sscanf( buf, " %f ", &x );
          j = this->GoToNextEOLData(j) +1;
          this->CellData[index]->InsertComponent( i-1, k, x); 
          }
        }
      }
    }
  // mccdo
  if ( index >= 0 )
    {
    if((this->CellData[index]->GetNumberOfTuples() == this->Mesh->GetNumberOfCells() )
      && (this->CellData[index]->GetNumberOfComponents() < 6) )
      {
        if (this->FirstArrayFlag == 0)
          {
          this->Mesh->GetCellData()->SetScalars(this->CellData[index]);
          }
        else
          {
          this->Mesh->GetCellData()->AddArray(this->CellData[index]);
          }
        this->FirstArrayFlag = 1;
        this->CellDataArraySelection->AddArray(this->CellData[ index ]->GetName());
        this->NumberOfCellFields++;
        this->CellData[index]->Delete();
        this->CellData[index] = 0;
      }
    }
  // mccdo
  return this->GoToNextSectionASCIIData(j);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataSinglePrecision(int ix)
{
  char buf[120];
  int index = -1; // mccdo 
  int j = ix + 1;
  j = this->GoToNextLeftParenData(j)+1;
  this->GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = this->GoToNextLeftParenData(j)+1;

  if ( this->DataPass == 1)
    {
    if ( this->IsCellZoneId(zid))
      {
      if ( this->IsNewVariable(ssid))
        {
        this->VariableIds->InsertValue(NumberOfVariables, ssid);
        this->VariableSizes->InsertValue(NumberOfVariables, size);
        this->NumberOfVariables++;
        }
      }
    }
  else
    {
    if ( this->IsCellZoneId(zid))
      {
      index = this->GetVariableIndex(ssid); // mccdo
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          this->CellData[index]->InsertComponent( i-1, k, 
            this->GetBinaryFloatData(j)); 
          j = j + 4;
          }
        }
      }
    }
  // mccdo
  if ( index >= 0 )
    {
    if((this->CellData[index]->GetNumberOfTuples() == this->Mesh->GetNumberOfCells() )
      && (this->CellData[index]->GetNumberOfComponents() < 6) )
      {
        if (this->FirstArrayFlag == 0)
          {
          this->Mesh->GetCellData()->SetScalars(this->CellData[index]);
          }
        else
          {
          this->Mesh->GetCellData()->AddArray(this->CellData[index]);
          }
        this->FirstArrayFlag = 1;
        this->CellDataArraySelection->AddArray(this->CellData[ index ]->GetName());
        this->NumberOfCellFields++;
        this->CellData[index]->Delete();
        this->CellData[index] = 0;
      }
    }
  // mccdo
  return this->GoToNextSectionSinglePrecisionData( j, "2300)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataDoublePrecision(int ix)
{
  char buf[120];
  int index = -1; //mccdo
  int j = ix + 1;
  j = this->GoToNextLeftParenData(j)+1;
  this->GetStringToNextRightParenData( j, buf );

  int ssid, zid, size, ntl, np, fi, li;
  sscanf( buf, " %d %d %d %d %d %d %d", 
    &ssid, &zid, &size, &ntl, &np, &fi, &li );
  j = this->GoToNextLeftParenData(j)+1;

  if ( this->DataPass == 1)
    {
    if ( this->IsCellZoneId(zid))
      {
      if ( this->IsNewVariable(ssid))
        {
        this->VariableIds->InsertValue(this->NumberOfVariables, ssid);
        this->VariableSizes->InsertValue(this->NumberOfVariables, size);
        this->NumberOfVariables++;
        }
      }
    }
  else
    {
    if ( this->IsCellZoneId(zid))
      {
      index = this->GetVariableIndex(ssid); // mccdo
      for (int i = fi; i <= li; i++)
        {
        for (int k = 0; k < size; k++)
          {
          this->CellData[index]->InsertComponent( i-1, k,
            this->GetBinaryDoubleData(j)); 
          j = j + 8;
          }
        }
      }
    }
  // mccdo
  if ( index >= 0 )
    {
    if((this->CellData[index]->GetNumberOfTuples() == this->Mesh->GetNumberOfCells() )
      && (this->CellData[index]->GetNumberOfComponents() < 6) )
      {
        if (this->FirstArrayFlag == 0)
          {
          this->Mesh->GetCellData()->SetScalars(this->CellData[index]);
          }
        else
          {
          this->Mesh->GetCellData()->AddArray(this->CellData[index]);
          }
        this->FirstArrayFlag = 1;
        this->CellDataArraySelection->AddArray(this->CellData[ index ]->GetName());
        this->NumberOfCellFields++;
        this->CellData[index]->Delete();
        this->CellData[index] = 0;
      }
    }
  // mccdo
  return this->GoToNextSectionSinglePrecisionData( j, "3300)");
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::SkipUnknownSinglePrecisionData(int ix, char buf[])
{
  return this->GoToNextSectionSinglePrecisionData( ix, buf);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::SkipUnknownDoublePrecisionData(int ix, char buf[])
{
  return this->GoToNextSectionDoublePrecisionData( ix, buf);
}


//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataUnknownASCII(int ix)
{
  int j = ix + 1;
  j = this->GoToNextLeftParenData(j)+1;
  j = this->GoToNextLeftParenData(j)+1;
  j = this->GoToNextRightParenData(j)+1;
  j = this->GoToNextRightParenData(j)+1;
  return j;
}

//-----------------------------------------------------------------------------
void vtkFLUENTReader::GetStringToNextRightParenData(int ix, char buf[] )
{
  // Copy contents between ( ) into buffer
  int j = ix;
  int k=0;
  while ( !(this->DataFileBuffer[j] == ')'))
    {
    buf[k] = this->DataFileBuffer[j];
    j++;
    k++;
    }
  buf[k] = 0;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::IsCellZoneId(int zi)
{
  int flag = 0;
  for (int i = 0; i < this->NumberOfCellZones; i++)
    {
    if ( zi == this->CellZones->GetValue(i))
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
  for (int i = 0; i < this->NumberOfVariables; i++)
    {
    if ( ssid == this->VariableIds->GetValue(i))
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
  for (int i = 0; i < this->NumberOfVariables; i++)
    {
    if ( ssid == this->VariableIds->GetValue(i))
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = this->DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = this->DataFileBuffer[ix+k];
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = this->DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = this->DataFileBuffer[ix+k];
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 7; k >= 0; k--)
      {
      mi.c[k] = this->DataFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 7; k++)
      {
      mi.c[7-k] = this->DataFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::IsASCIICharacterHexDigit(int ix)
{
  if ( (this->CaseFileBuffer[ix] >= 0x30) 
    && (this->CaseFileBuffer[ix] <= 0x39))
    {
    return 1;
    }
  else if ( (this->CaseFileBuffer[ix] >= 0x41) 
    && (this->CaseFileBuffer[ix] <= 0x46))
    {
    return 1;
    }
  else if ( (this->CaseFileBuffer[ix] >= 0x61) 
    && (this->CaseFileBuffer[ix] <= 0x66))
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
  while (! this->IsASCIICharacterHexDigit(i))
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextRightParen(int ix)
{
  int i = ix;
  while (this->CaseFileBuffer[i] != ')' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextLeftParen(int ix)
{
  int i = ix;
  while (this->CaseFileBuffer[i] != '(' )
    {
    i++;
    }
  return i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextEOL(int ix)
{
  int i = ix;
  while (this->CaseFileBuffer[i] != 0x0a )
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
  while ( !((level == 0) && (this->CaseFileBuffer[i] == ')')))
    {
    if (this->CaseFileBuffer[i] == '(')
      {
      level++;
      }
    if (this->CaseFileBuffer[i] == ')')
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
  while ( !(this->CaseFileBuffer[j] == ')'))
    {
    buf[k] = this->CaseFileBuffer[j];
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
  while ( !((this->CaseFileBuffer[j] == 0x0a)
    ||(this->CaseFileBuffer[j] == ')')))
    {
    buf[k] = this->CaseFileBuffer[j];
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
    j = this->GoToNextASCIIHexDigit(j);
    //cout << "i = " << i << ", et = " << this->CaseFileBuffer[j] << endl;
    c[0] = this->CaseFileBuffer[j];
    c[1] = 0;
    this->CellTypes->InsertValue(i, atoi(c));
    j++;
    }
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionSinglePrecision(int ix, char buf[])
{
  int i = ix + 1;
  while( !((this->CaseFileBuffer[i] == buf[0]) 
    && (this->CaseFileBuffer[i+1] == buf[1]) 
    && (this->CaseFileBuffer[i+2] == buf[2]) 
    && (this->CaseFileBuffer[i+3] == buf[3]) 
    && (this->CaseFileBuffer[i+4] == buf[4]) ))
    {
    i++;
    }
  return i+4;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextSectionDoublePrecision(int ix, char buf[])
{
  int i = ix + 1;
  while ( !((this->CaseFileBuffer[i] == buf[0]) 
    && (this->CaseFileBuffer[i+1] == buf[1]) 
    && (this->CaseFileBuffer[i+2] == buf[2]) 
    && (this->CaseFileBuffer[i+3] == buf[3]) 
    && (this->CaseFileBuffer[i+4] == buf[4])))
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = this->CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = this->CaseFileBuffer[ix+k];
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 3; k >= 0; k--)
      {
      mi.c[k] = this->CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 3; k++)
      {
      mi.c[3-k] = this->CaseFileBuffer[ix+k];
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

  if ( this->LittleEndianFlag == 1)
    {
    for (int k = 7; k >= 0; k--)
      {
      mi.c[k] = this->CaseFileBuffer[ix+k];
      }
    }
  else
    {
    for (int k = 0; k <= 7; k++)
      {
      mi.c[7-k] = this->CaseFileBuffer[ix+k];
      }
    }
  return mi.i;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GetAsciiInteger(int ix)
{
  int j = ix;
  int first = this->GoToNextASCIIHexDigit(j);
  j = first;
  while ( this->IsASCIICharacterHexDigit(++j));
  int second = j-1;
  int nchar = second-first+1;
  char *buf;
  buf = new char[nchar+1];
  buf[nchar] = 0;
  j = first;

  for (int i = 0; i < nchar; i++)
    {
    buf[i] = this->CaseFileBuffer[j];
    j++;
    }
  return atoi(buf);
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoPastAsciiInteger(int ix)
{
  int j = ix;
  int first = this->GoToNextASCIIHexDigit(j);
  j = first;
  while(this->IsASCIICharacterHexDigit(++j));
  return j;
}

//-----------------------------------------------------------------------------
int vtkFLUENTReader::GoToNextEOLData(int ix)
{
  int i = ix;
  while(this->DataFileBuffer[i] != 0x0a )
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
  while ( !((this->DataFileBuffer[j] == 0x0a)
    ||(this->DataFileBuffer[j] == ')'))) 
    {
    buf[k] = this->DataFileBuffer[j];
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
  // mccdo this->CellFaces = vtkIntArray::New();
  // mccdo this->CellFacesClean = vtkIntArray::New();
  // mccdo this->CellFacesClean->SetNumberOfComponents(6);
  this->FaceTypes = vtkIntArray::New();
  this->FaceNodes = vtkIntArray::New();
  this->FaceNodes->SetNumberOfComponents(4);
  this->FaceCells = vtkIntArray::New();
  this->FaceCells->SetNumberOfComponents(2);
  // mccdo this->FaceParents = vtkIntArray::New();
  // mccdo this->FaceParents->SetNumberOfComponents(2);
  // mccdo this->PeriodicShadowFaces = vtkIntArray::New();
  // mccdo this->PeriodicShadowFaces->SetNumberOfComponents(2);
  this->FaceTreesNumberOfKids = vtkIntArray::New();
  this->FaceTreesKids = vtkIntArray::New();
  this->FaceTreesKidsIndex = vtkIntArray::New();
  // mccdo this->CellTreesNumberOfKids = vtkIntArray::New();
  // mccdo this->CellTreesKids = vtkIntArray::New();
  // mccdo this->CellTreesKidsIndex = vtkIntArray::New();
  // mccdo this->FaceTreeParentFaceId0 = vtkIntArray::New();
  // mccdo this->FaceTreeParentFaceId1 = vtkIntArray::New();
  this->FaceTreeParentTable = vtkIntArray::New();
  // mccdo this->CellTreeParentCellId0 = vtkIntArray::New();
  // mccdo this->CellTreeParentCellId1 = vtkIntArray::New();
  // mccdo this->NCGFaceChild = vtkIntArray::New();
  // mccdo this->NCGFaceParent = vtkIntArray::New();
  // mccdo this->CellNumberOfFaces = vtkIntArray::New();
  // mccdo this->FaceParentFlags = vtkIntArray::New();
  // mccdo this->CellIndex = vtkIntArray::New();
  // mccdo this->InterfaceFaceChildFlags = vtkIntArray::New();
  // mccdo this->FaceParentsChildren = vtkIntArray::New();
  // mccdo this->NCGFaceChildFlags = vtkIntArray::New();
  // mccdo this->CellParentFlags = vtkIntArray::New();
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

  this->ObjectsFlag = 1;
  this->FirstArrayFlag = 0;

}


//-----------------------------------------------------------------------------
void vtkFLUENTReader::DeleteVTKObjects(void)
{
  // mccdo delete [] CaseFileBuffer;
  // mccdo delete [] DataFileBuffer;
  // mccdo this->Points->Delete();
  // mccdo this->CellTypes->Delete();
  // mccdo this->CellFaces->Delete();
  // mccdo this->CellFacesClean->Delete();

  // mccdo this->FaceTypes->Delete();
  // mccdo this->FaceNodes->Delete();
  // mccdo this->FaceCells->Delete();
  // mccdo this->FaceParents->Delete();
  // mccdo this->PeriodicShadowFaces->Delete();
  // mccdo this->FaceTreesNumberOfKids->Delete();
  // mccdo this->FaceTreesKids->Delete();
  // mccdo this->FaceTreesKidsIndex->Delete();
  // mccdo this->CellTreesNumberOfKids->Delete();
  // mccdo this->CellTreesKids->Delete();
  // mccdo this->CellTreesKidsIndex->Delete();
  // mccdo this->FaceTreeParentTable->Delete();
  // mccdo this->CellTreeParentCellId0->Delete();
  // mccdo this->CellTreeParentCellId1->Delete();

  // mccdo this->NCGFaceParent->Delete();
  // mccdo this->CellNumberOfFaces->Delete();
  // mccdo this->FaceParentFlags->Delete();
  // mccdo this->CellIndex->Delete();
  // mccdo this->InterfaceFaceChildFlags->Delete();
  // mccdo this->FaceParentsChildren->Delete();
  // mccdo this->NCGFaceChildFlags->Delete();
  // mccdo this->CellParentFlags->Delete();
  //this->CellDataArraySelection->Delete(); // mccdo
  this->ATriangle->Delete();
  this->AQuad->Delete();
  this->ATetra->Delete();
  this->APyramid->Delete();
  this->AWedge->Delete();
  this->AHexahedron->Delete();
  this->CellZones->Delete();
  this->VariableIds->Delete();
  this->VariableSizes->Delete();
  this->ObjectsFlag = 0;
}
