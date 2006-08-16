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
#include "vtkMultiBlockDataSet.h"
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
#include "vtkConvexPointSet.h"
#include <fstream>
#include <sstream>
using vtkstd::istringstream;

vtkCxxRevisionMacro(vtkFLUENTReader, "$Revision$");
vtkStandardNewMacro(vtkFLUENTReader);

//----------------------------------------------------------------------------
vtkFLUENTReader::vtkFLUENTReader()
{cout<<"Constructor"<<endl;
  this->SetNumberOfInputPorts(0);
  this->FileName  = NULL;
  this->Points = vtkPoints::New();
  this->Triangle = vtkTriangle::New();
  this->Tetra = vtkTetra::New();
  this->Quad = vtkQuad::New();
  this->Hexahedron = vtkHexahedron::New();
  this->Pyramid = vtkPyramid::New();
  this->Wedge = vtkWedge::New();
  this->ConvexPointSet = vtkConvexPointSet::New();
cout<<"ConstructorEnd"<<endl;
}

//----------------------------------------------------------------------------
vtkFLUENTReader::~vtkFLUENTReader()
{cout<<"deConstructor"<<endl;
  Points->Delete();
  Triangle->Delete();
  Tetra->Delete();
  Quad->Delete();
  Hexahedron->Delete();
  Pyramid->Delete();
  Wedge->Delete();
  ConvexPointSet->Delete();
cout<<"deConstructorEnd"<<endl;
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{cout<<"RequestData"<<endl;
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkMultiBlockDataSet::COMPOSITE_DATA_SET()));

  output->SetNumberOfDataSets(0, CellZones.size());

  vtkUnstructuredGrid *Grid[CellZones.size()];

  for(int test=0; test < CellZones.size(); test++)
    Grid[test] = vtkUnstructuredGrid::New();

  for (int i = 0; i < Cells.size(); i++)
    {
    int location = distance(CellZones.begin(),find(CellZones.begin(), CellZones.end(), Cells[i].zone));

    if (Cells[i].type == 1 )
      {
      for (int j = 0; j < 3; j++)
        {
        Triangle->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Triangle->GetCellType(),Triangle->GetPointIds());
      }
    else if (Cells[i].type == 2 )
      {
      for (int j = 0; j < 4; j++)
        {
        Tetra->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Tetra->GetCellType(),Tetra->GetPointIds());
      }
    else if (Cells[i].type == 3 )
      {
      for (int j = 0; j < 4; j++)
        {
        Quad->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Quad->GetCellType(),Quad->GetPointIds());
      }
    else if (Cells[i].type == 4 )
      {
      for (int j = 0; j < 8; j++)
        {
        Hexahedron->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Hexahedron->GetCellType(),Hexahedron->GetPointIds());
      }
    else if (Cells[i].type == 5 )
      {
      for (int j = 0; j < 5; j++)
        {
        Pyramid->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Pyramid->GetCellType(),Pyramid->GetPointIds());
      }
    else if (Cells[i].type == 6 )
      {
      for (int j = 0; j < 6; j++)
        {
        Wedge->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(Wedge->GetCellType(),Wedge->GetPointIds());
      }
    else if (Cells[i].type == 7 )
      {
      ConvexPointSet->GetPointIds()->SetNumberOfIds(Cells[i].nodes.size());
      for (int j = 0; j < Cells[i].nodes.size(); j++)
        {
        ConvexPointSet->GetPointIds()->SetId( j, Cells[i].nodes[j]);
        }
      Grid[location]->InsertNextCell(ConvexPointSet->GetCellType(),ConvexPointSet->GetPointIds());
      }
    }
  Cells.clear();

  //Scalar Data
  for (int l = 0; l < ScalarDataChunks.size(); l++)
    {
    int location = distance(CellZones.begin(),find(CellZones.begin(), CellZones.end(), ScalarDataChunks[l].zoneId));
    vtkDoubleArray *v = vtkDoubleArray::New();
    for (int m = 0; m < ScalarDataChunks[l].scalarData.size(); m++)
      {
      v->InsertValue(m, ScalarDataChunks[l].scalarData[m]);
      }
    v->SetName(ScalarVariableNames[l/CellZones.size()].c_str());
    Grid[location]->GetCellData()->AddArray(v);
    v->Delete();
    }

  //Vector Data
  for (int l = 0; l < VectorDataChunks.size(); l++)
    {
    int location = distance(CellZones.begin(),find(CellZones.begin(), CellZones.end(), VectorDataChunks[l].zoneId));
    vtkDoubleArray *v = vtkDoubleArray::New();
    v->SetNumberOfComponents(3);
    for (int m = 0; m < VectorDataChunks[l].iComponentData.size(); m++)
      {
      v->InsertComponent(m, 0, VectorDataChunks[l].iComponentData[m]);
      v->InsertComponent(m, 1, VectorDataChunks[l].jComponentData[m]);
      v->InsertComponent(m, 2, VectorDataChunks[l].kComponentData[m]);
      }
    v->SetName(VectorVariableNames[l/CellZones.size()].c_str());
    Grid[location]->GetCellData()->AddArray(v);
    v->Delete();
    }

  for(int addTo = 0; addTo < CellZones.size(); addTo++)
    {
    Grid[addTo]->SetPoints(Points);
    output->SetDataSet(0, addTo, Grid[addTo]);
    Grid[addTo]->Delete();
    }
  cout<<"RequestDataEnd"<<endl;
  return 1;
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PrintSelf(ostream& os, vtkIndent indent)
{cout<<"PrintSelf"<<endl;
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << endl;

  os << indent << "Number Of Cells: " << this->NumberOfCells << endl;
cout<<"PrintSelfEnd"<<endl;
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{cout<<"RequestInfo"<<endl;
  OpenCaseFile(this->FileName);
  OpenDataFile(this->FileName);
cout<<"Case Opened"<<endl;
  ParseCaseFile();  // Reads Necessary Information from the .cas file.

cout<<"Case Parsed"<<endl;
  CleanCells();  //  Removes unnecessary faces from the cells.

cout<<"Cells Cleaned"<<endl;
  PopulateCellNodes();

cout<<"Cells populated"<<endl;
  LoadVariableNames();
cout<<"Names Loaded"<<endl;

  GetNumberOfCellZones();
cout<<CellZones.size()<<endl;
  NumberOfScalars = 0;
  NumberOfVectors = 0;
  ParseDataFile();
cout<<SubSectionIds.size()<<endl;
  this->CellDataArraySelection = vtkDataArraySelection::New();
  for (int i = 0; i < SubSectionIds.size(); i++)
    {
cout<<i<<" : "<<SubSectionIds[i]<<" : "<<SubSectionSize[i]<<" : "<<VariableNames[SubSectionIds[i]]<<endl;
    if (SubSectionSize[i] == 1)
      {
      cout<<"Scalar"<<endl;
      this->CellDataArraySelection->AddArray(VariableNames[SubSectionIds[i]].c_str());
      ScalarVariableNames.push_back(VariableNames[SubSectionIds[i]]);
      ScalarSubSectionIds.push_back(SubSectionIds[i]);
      }
    else if (SubSectionSize[i] == 3)
      {
      cout<<"Vector"<<endl;
      this->CellDataArraySelection->AddArray(VariableNames[SubSectionIds[i]].c_str());
      VectorVariableNames.push_back(VariableNames[SubSectionIds[i]]);
      VectorSubSectionIds.push_back(SubSectionIds[i]);
      }
 cout<<this->CellDataArraySelection->GetArrayName(i)<<endl;
    }
cout<<"RequestInfoEnd"<<endl;
return 1;
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::OpenCaseFile(const char *filename)
{
  this->FluentCaseFile.open(filename);

  if (this->FluentCaseFile.is_open())
    {
    cout << "Successfully opened " << filename << endl;
    return 1;
    }
  else
    {
    cout << "Could not open " << filename << endl;
    return 0;
    }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
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

//----------------------------------------------------------------------------
int vtkFLUENTReader::OpenDataFile(const char *filename)
{
  vtkstd::string dfilename(filename);
  dfilename.erase(dfilename.length()-3, 3);
  dfilename.append("dat");

  FluentDataFile.open(dfilename.c_str());

  if (FluentDataFile.is_open())
    {
    cout << "Successfully opened " << dfilename << endl;
    return 1;
    }
  else
    {
    cout << "Could not open " << dfilename << endl;
    return 0;
    }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetCaseChunk ()
{
cout<<"Get Case Chunk"<<endl;
  CaseBuffer.clear();  // Clear buffer

  //
  // Look for beginning of chunk
  //
  while(FluentCaseFile.peek() != '(')
    {
    char c = FluentCaseFile.get();
    if (FluentCaseFile.eof())
      {
      return 0;
      }
    }

  //
  // Figure out whether this is a binary or ascii chunk.
  // If the index is 3 digits or more, then binary, otherwise ascii.
  //
  vtkstd::string index;
  while(FluentCaseFile.peek() != ' ')
    {
    index.push_back(FluentCaseFile.peek());
    CaseBuffer.push_back(FluentCaseFile.get());
    if (FluentCaseFile.eof())
      {
      return 0;
      }
    }

  index.erase(0,1);  // Get rid of the "("

  //
  //  Grab the chunk and put it in buffer.
  //  You have to look for the end of section vtkstd::string if it is
  //  a binary chunk.
  //

  if (index.size() > 2)
    {  // Binary Chunk
    char end[120];
    strcpy(end, "End of Binary Section   ");
    strcat(end, index.c_str());
    strcat(end, ")");

    // Load the case buffer enough to start comparing to the end vtkstd::string.
    while (CaseBuffer.size() < strlen(end))
      {
      CaseBuffer.push_back(FluentCaseFile.get());
      }

    while ( CaseBuffer.compare(CaseBuffer.size()-strlen(end), strlen(end), end) )
      {
      CaseBuffer.push_back(FluentCaseFile.get());
      }

    }
  else
    {  // Ascii Chunk
    int level = 0;
    while ((FluentCaseFile.peek() != ')') || (level != 0) )
      {
      CaseBuffer.push_back(FluentCaseFile.get());
      if (CaseBuffer.at(CaseBuffer.length()-1) == '(')
        {
        level++;
        }
      if (CaseBuffer.at(CaseBuffer.length()-1) == ')')
        {
        level--;
        }
      if (FluentCaseFile.eof())
        {
        return 0;
        }
      }
    CaseBuffer.push_back(FluentCaseFile.get());
    }
  return 1;
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetCaseIndex()
{
  vtkstd::string sindex;

  int i = 1;
  while (CaseBuffer.at(i) != ' ')
    {
    sindex.push_back(CaseBuffer.at(i++));
    }
  return atoi(sindex.c_str());
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNumberOfCellZones()
{
  int match;

  for (int i = 0; i < Cells.size(); i++)
    {
    if (CellZones.size() == 0)
      {
      CellZones.push_back(Cells[i].zone);
      }
    else 
      {
      match = 0;
      for (int j = 0; j < CellZones.size(); j++)
        {
        if (CellZones[j] == Cells[i].zone)
          {
          match = 1;
          }
        }
        if (match == 0)
          {
          CellZones.push_back(Cells[i].zone);
          }
      }
    }
}


//----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataIndex()
{
  vtkstd::string sindex;

  int i = 1;
  while (DataBuffer.at(i) != ' ')
    {
    sindex.push_back(DataBuffer.at(i++));
    }
  return atoi(sindex.c_str());
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataChunk ()
{
  DataBuffer.clear();  // Clear buffer
  //
  // Look for beginning of chunk
  //
cout<<FluentDataFile.peek()<<endl;
  while(FluentDataFile.peek() != '(')
    {
    char c = FluentDataFile.get();
cout<<"gdc:"<<c<<endl;
    if (FluentDataFile.eof())
      {
      return 0;
      }
    }

  //
  // Figure out whether this is a binary or ascii chunk.
  // If the index is 3 digits or more, then binary, otherwise ascii.
  //
  vtkstd::string index;
  while(FluentDataFile.peek() != ' ')
    {
    index.push_back(FluentDataFile.peek());
    DataBuffer.push_back(FluentDataFile.get());
    if (FluentDataFile.eof())
      {
      return 0;
      }
    }

  index.erase(0,1);  // Get rid of the "("

  //
  //  Grab the chunk and put it in buffer.
  //  You have to look for the end of section vtkstd::string if it is
  //  a binary chunk.
  //
  if (index.size() > 3)
    {  // Binary Chunk
    char end[120];
    strcpy(end, "End of Binary Section   ");
    strcat(end, index.c_str());
    strcat(end, ")");

    // Load the data buffer enough to start comparing to the end vtkstd::string.
    while (DataBuffer.size() < strlen(end))
      {
      DataBuffer.push_back(FluentDataFile.get());
      }

    while ( DataBuffer.compare(DataBuffer.size()-strlen(end), strlen(end), end) )
      {
      DataBuffer.push_back(FluentDataFile.get());
      }

    }
  else
    {  // Ascii Chunk
    int level = 0;
    while ((FluentDataFile.peek() != ')') || (level != 0) )
      {
      DataBuffer.push_back(FluentDataFile.get());
      if (DataBuffer.at(DataBuffer.length()-1) == '(')
        {
        level++;
        }
      if (DataBuffer.at(DataBuffer.length()-1) == ')')
        {
        level--;
        }
      if (FluentDataFile.eof())
        {
        return 0;
        }
      }
    DataBuffer.push_back(FluentDataFile.get());
    }

  return 1;
}


void vtkFLUENTReader::LoadVariableNames()
{
  VariableNames[1]  = "PRESSURE";
  VariableNames[2]  = "MOMENTUM";
  VariableNames[3]  = "TEMPERATURE";
  VariableNames[4]  = "ENTHALPY";
  VariableNames[5]  = "TKE";
  VariableNames[6]  = "TED";
  VariableNames[7]  = "SPECIES";
  VariableNames[8]  = "G";
  VariableNames[9]  = "WSWIRL";
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
  VariableNames[56] = "DO_IWX=56";
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
  VariableNames[182] = "POLLUT7";
  VariableNames[183] = "POLLUT8";
  VariableNames[184] = "POLLUT9";
  VariableNames[185] = "POLLUT10";
  VariableNames[186] = "POLLUT11";
  VariableNames[187] = "POLLUT12";
  VariableNames[188] = "POLLUT13";
  VariableNames[190] = "SV_T_AUX";
  VariableNames[191] = "SV_T_AP_AUX";
  VariableNames[192] = "TOTAL_PRESSURE";
  VariableNames[193] = "TOTAL_TEMPERATURE";
  VariableNames[194] = "NRBC_DC";
  VariableNames[195] = "DP_TMFR";
  
  
  VariableNames[200] = "Y_00"; 
  VariableNames[201] = "Y_01"; 
  VariableNames[202] = "Y_02"; 
  VariableNames[203] = "Y_03"; 
  VariableNames[204] = "Y_04"; 
  VariableNames[205] = "Y_05"; 
  VariableNames[206] = "Y_06"; 
  VariableNames[207] = "Y_07"; 
  VariableNames[208] = "Y_08"; 
  VariableNames[209] = "Y_09"; 
  VariableNames[210] = "Y_10"; 
  VariableNames[211] = "Y_11"; 
  VariableNames[212] = "Y_12"; 
  VariableNames[213] = "Y_13"; 
  VariableNames[214] = "Y_14"; 
  VariableNames[215] = "Y_15"; 
  VariableNames[216] = "Y_16"; 
  VariableNames[217] = "Y_17"; 
  VariableNames[218] = "Y_18"; 
  VariableNames[219] = "Y_19"; 
  VariableNames[220] = "Y_20"; 
  VariableNames[221] = "Y_21"; 
  VariableNames[222] = "Y_22"; 
  VariableNames[223] = "Y_23"; 
  VariableNames[224] = "Y_24"; 
  VariableNames[225] = "Y_25"; 
  VariableNames[226] = "Y_26"; 
  VariableNames[227] = "Y_27"; 
  VariableNames[228] = "Y_28"; 
  VariableNames[229] = "Y_29"; 
  VariableNames[230] = "Y_30"; 
  VariableNames[231] = "Y_31"; 
  VariableNames[232] = "Y_32"; 
  VariableNames[233] = "Y_33"; 
  VariableNames[234] = "Y_34"; 
  VariableNames[235] = "Y_35"; 
  VariableNames[236] = "Y_36"; 
  VariableNames[237] = "Y_37"; 
  VariableNames[238] = "Y_38"; 
  VariableNames[239] = "Y_39"; 
  VariableNames[240] = "Y_40"; 
  VariableNames[241] = "Y_41"; 
  VariableNames[242] = "Y_42"; 
  VariableNames[243] = "Y_43"; 
  VariableNames[244] = "Y_44"; 
  VariableNames[245] = "Y_45"; 
  VariableNames[246] = "Y_46"; 
  VariableNames[247] = "Y_47"; 
  VariableNames[248] = "Y_48"; 
  VariableNames[249] = "Y_49"; 

  VariableNames[250] = "Y_M1_00"; 
  VariableNames[251] = "Y_M1_01"; 
  VariableNames[252] = "Y_M1_02"; 
  VariableNames[253] = "Y_M1_03"; 
  VariableNames[254] = "Y_M1_04"; 
  VariableNames[255] = "Y_M1_05"; 
  VariableNames[256] = "Y_M1_06"; 
  VariableNames[257] = "Y_M1_07"; 
  VariableNames[258] = "Y_M1_08"; 
  VariableNames[259] = "Y_M1_09"; 
  VariableNames[260] = "Y_M1_10"; 
  VariableNames[261] = "Y_M1_11"; 
  VariableNames[262] = "Y_M1_12"; 
  VariableNames[263] = "Y_M1_13"; 
  VariableNames[264] = "Y_M1_14"; 
  VariableNames[265] = "Y_M1_15"; 
  VariableNames[266] = "Y_M1_16"; 
  VariableNames[267] = "Y_M1_17"; 
  VariableNames[268] = "Y_M1_18"; 
  VariableNames[269] = "Y_M1_19"; 
  VariableNames[270] = "Y_M1_20"; 
  VariableNames[271] = "Y_M1_21"; 
  VariableNames[272] = "Y_M1_22"; 
  VariableNames[273] = "Y_M1_23"; 
  VariableNames[274] = "Y_M1_24"; 
  VariableNames[275] = "Y_M1_25"; 
  VariableNames[276] = "Y_M1_26"; 
  VariableNames[277] = "Y_M1_27"; 
  VariableNames[278] = "Y_M1_28"; 
  VariableNames[279] = "Y_M1_29"; 
  VariableNames[280] = "Y_M1_30"; 
  VariableNames[281] = "Y_M1_31"; 
  VariableNames[282] = "Y_M1_32"; 
  VariableNames[283] = "Y_M1_33"; 
  VariableNames[284] = "Y_M1_34"; 
  VariableNames[285] = "Y_M1_35"; 
  VariableNames[286] = "Y_M1_36"; 
  VariableNames[287] = "Y_M1_37"; 
  VariableNames[288] = "Y_M1_38"; 
  VariableNames[289] = "Y_M1_39"; 
  VariableNames[290] = "Y_M1_40"; 
  VariableNames[291] = "Y_M1_41"; 
  VariableNames[292] = "Y_M1_42"; 
  VariableNames[293] = "Y_M1_43"; 
  VariableNames[294] = "Y_M1_44"; 
  VariableNames[295] = "Y_M1_45"; 
  VariableNames[296] = "Y_M1_46"; 
  VariableNames[297] = "Y_M1_47"; 
  VariableNames[298] = "Y_M1_48"; 
  VariableNames[299] = "Y_M1_49"; 

  VariableNames[300] = "Y_M2_00"; 
  VariableNames[301] = "Y_M2_01"; 
  VariableNames[302] = "Y_M2_02"; 
  VariableNames[303] = "Y_M2_03"; 
  VariableNames[304] = "Y_M2_04"; 
  VariableNames[305] = "Y_M2_05"; 
  VariableNames[306] = "Y_M2_06"; 
  VariableNames[307] = "Y_M2_07"; 
  VariableNames[308] = "Y_M2_08"; 
  VariableNames[309] = "Y_M2_09"; 
  VariableNames[310] = "Y_M2_10"; 
  VariableNames[311] = "Y_M2_11"; 
  VariableNames[312] = "Y_M2_12"; 
  VariableNames[313] = "Y_M2_13"; 
  VariableNames[314] = "Y_M2_14"; 
  VariableNames[315] = "Y_M2_15"; 
  VariableNames[316] = "Y_M2_16"; 
  VariableNames[317] = "Y_M2_17"; 
  VariableNames[318] = "Y_M2_18"; 
  VariableNames[319] = "Y_M2_19"; 
  VariableNames[320] = "Y_M2_20"; 
  VariableNames[321] = "Y_M2_21"; 
  VariableNames[322] = "Y_M2_22"; 
  VariableNames[323] = "Y_M2_23"; 
  VariableNames[324] = "Y_M2_24"; 
  VariableNames[325] = "Y_M2_25"; 
  VariableNames[326] = "Y_M2_26"; 
  VariableNames[327] = "Y_M2_27"; 
  VariableNames[328] = "Y_M2_28"; 
  VariableNames[329] = "Y_M2_29"; 
  VariableNames[330] = "Y_M2_30"; 
  VariableNames[331] = "Y_M2_31"; 
  VariableNames[332] = "Y_M2_32"; 
  VariableNames[333] = "Y_M2_33"; 
  VariableNames[334] = "Y_M2_34"; 
  VariableNames[335] = "Y_M2_35"; 
  VariableNames[336] = "Y_M2_36"; 
  VariableNames[337] = "Y_M2_37"; 
  VariableNames[338] = "Y_M2_38"; 
  VariableNames[339] = "Y_M2_39"; 
  VariableNames[340] = "Y_M2_40"; 
  VariableNames[341] = "Y_M2_41"; 
  VariableNames[342] = "Y_M2_42"; 
  VariableNames[343] = "Y_M2_43"; 
  VariableNames[344] = "Y_M2_44"; 
  VariableNames[345] = "Y_M2_45"; 
  VariableNames[346] = "Y_M2_46"; 
  VariableNames[347] = "Y_M2_47"; 
  VariableNames[348] = "Y_M2_48"; 
  VariableNames[349] = "Y_M2_49"; 

  VariableNames[350] = "DR_SURF_00"; 
  VariableNames[351] = "DR_SURF_01"; 
  VariableNames[352] = "DR_SURF_02"; 
  VariableNames[353] = "DR_SURF_03"; 
  VariableNames[354] = "DR_SURF_04"; 
  VariableNames[355] = "DR_SURF_05"; 
  VariableNames[356] = "DR_SURF_06"; 
  VariableNames[357] = "DR_SURF_07"; 
  VariableNames[358] = "DR_SURF_08"; 
  VariableNames[359] = "DR_SURF_09"; 
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

  VariableNames[450] = "DPMS_Y_00"; 
  VariableNames[451] = "DPMS_Y_01"; 
  VariableNames[452] = "DPMS_Y_02"; 
  VariableNames[453] = "DPMS_Y_03"; 
  VariableNames[454] = "DPMS_Y_04"; 
  VariableNames[455] = "DPMS_Y_05"; 
  VariableNames[456] = "DPMS_Y_06"; 
  VariableNames[457] = "DPMS_Y_07"; 
  VariableNames[458] = "DPMS_Y_08"; 
  VariableNames[459] = "DPMS_Y_09"; 
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
  VariableNames[527] = "LKE";
  VariableNames[528] = "LKE_M1";
  VariableNames[529] = "LKE_M2";
  VariableNames[530] = "SHELL_CELL_T";
  VariableNames[531] = "SHELL_FACE_T";
  VariableNames[532] = "SHELL_CELL_ENERGY_M1";
  VariableNames[533] = "SHELL_CELL_ENERGY_M2";
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
  VariableNames[580] = "MU_TURB_L";
  VariableNames[581] = "MU_TURB_S";
  VariableNames[582] = "TKE_TRANS";
  VariableNames[583] = "TKE_TRANS_M1";
  VariableNames[584] = "TKE_TRANS_M2";
  VariableNames[585] = "MU_TURB_W";
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
  VariableNames[620] = "MICRO_MIX_FVAR1 "; 
  VariableNames[621] = "MICRO_MIX_FVAR2 "; 
  VariableNames[622] = "MICRO_MIX_FVAR3 "; 
  VariableNames[623] = "MICRO_MIX_FVAR1_M1 "; 
  VariableNames[624] = "MICRO_MIX_FVAR2_M1 "; 
  VariableNames[625] = "MICRO_MIX_FVAR3_M1 "; 
  VariableNames[626] = "MICRO_MIX_FVAR1_M2 "; 
  VariableNames[627] = "MICRO_MIX_FVAR2_M2 "; 
  VariableNames[628] = "MICRO_MIX_FVAR3_M2 "; 
  VariableNames[630] = "SCAD_LES "; 
  VariableNames[635] = "UFLA_Y    "; 
  VariableNames[636] = "UFLA_Y_M1 "; 
  VariableNames[637] = "UFLA_Y_M2 "; 
  VariableNames[645] = "CREV_MASS";
  VariableNames[646] = "CREV_ENRG";
  VariableNames[647] = "CREV_MOM";
  VariableNames[650] = "ACOUSTICS_MODEL";
  VariableNames[651] = "AC_RECEIVERS_DATA";
  VariableNames[652] = "SV_DPDT_RMS"; 
  VariableNames[653] = "SV_PRESSURE_M1"; 
  VariableNames[654] = "AC_PERIODIC_INDEX"; 
  VariableNames[655] = "AC_PERIODIC_PS";
  VariableNames[656] = "AC_F_NORMAL";
  VariableNames[657] = "AC_F_CENTROID";
  VariableNames[660] = "IGNITE";
  VariableNames[661] = "IGNITE_M1";
  VariableNames[662] = "IGNITE_M2";
  VariableNames[663] = "IGNITE_RATE";

  VariableNames[680] = "WALL_SHEAR_MEAN";
  VariableNames[681] = "UV_MEAN";
  VariableNames[682] = "UW_MEAN";
  VariableNames[683] = "VW_MEAN";
  VariableNames[684] = "UT_MEAN";
  VariableNames[685] = "VT_MEAN";
  VariableNames[686] = "WT_MEAN";
  VariableNames[687] = "BOUNDARY_HEAT_FLUX_MEAN";

  VariableNames[700] = "UDS_00"; 
  VariableNames[701] = "UDS_01"; 
  VariableNames[702] = "UDS_02"; 
  VariableNames[703] = "UDS_03"; 
  VariableNames[704] = "UDS_04"; 
  VariableNames[705] = "UDS_05"; 
  VariableNames[706] = "UDS_06"; 
  VariableNames[707] = "UDS_07"; 
  VariableNames[708] = "UDS_08"; 
  VariableNames[709] = "UDS_09"; 
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

  VariableNames[750] = "UDS_M1_00"; 
  VariableNames[751] = "UDS_M1_01"; 
  VariableNames[752] = "UDS_M1_02"; 
  VariableNames[753] = "UDS_M1_03"; 
  VariableNames[754] = "UDS_M1_04"; 
  VariableNames[755] = "UDS_M1_05"; 
  VariableNames[756] = "UDS_M1_06"; 
  VariableNames[757] = "UDS_M1_07"; 
  VariableNames[758] = "UDS_M1_08"; 
  VariableNames[759] = "UDS_M1_09"; 
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

  VariableNames[800] = "UDS_M2_00"; 
  VariableNames[801] = "UDS_M2_01"; 
  VariableNames[802] = "UDS_M2_02"; 
  VariableNames[803] = "UDS_M2_03"; 
  VariableNames[804] = "UDS_M2_04"; 
  VariableNames[805] = "UDS_M2_05"; 
  VariableNames[806] = "UDS_M2_06"; 
  VariableNames[807] = "UDS_M2_07"; 
  VariableNames[808] = "UDS_M2_08"; 
  VariableNames[809] = "UDS_M2_09"; 
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

  VariableNames[850] = "DPMS_DS_Y_00"; 
  VariableNames[851] = "DPMS_DS_Y_01"; 
  VariableNames[852] = "DPMS_DS_Y_02"; 
  VariableNames[853] = "DPMS_DS_Y_03"; 
  VariableNames[854] = "DPMS_DS_Y_04"; 
  VariableNames[855] = "DPMS_DS_Y_05"; 
  VariableNames[856] = "DPMS_DS_Y_06"; 
  VariableNames[857] = "DPMS_DS_Y_07"; 
  VariableNames[858] = "DPMS_DS_Y_08"; 
  VariableNames[859] = "DPMS_DS_Y_09"; 
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

  VariableNames[920] = "DPMS_DS_SURFACE_SPECIES_00"; 
  VariableNames[921] = "DPMS_DS_SURFACE_SPECIES_01"; 
  VariableNames[922] = "DPMS_DS_SURFACE_SPECIES_02"; 
  VariableNames[923] = "DPMS_DS_SURFACE_SPECIES_03"; 
  VariableNames[924] = "DPMS_DS_SURFACE_SPECIES_04"; 
  VariableNames[925] = "DPMS_DS_SURFACE_SPECIES_05"; 
  VariableNames[926] = "DPMS_DS_SURFACE_SPECIES_06"; 
  VariableNames[927] = "DPMS_DS_SURFACE_SPECIES_07"; 
  VariableNames[928] = "DPMS_DS_SURFACE_SPECIES_08"; 
  VariableNames[929] = "DPMS_DS_SURFACE_SPECIES_09"; 
  VariableNames[930] = "DPMS_DS_SURFACE_SPECIES_10"; 
  VariableNames[931] = "DPMS_DS_SURFACE_SPECIES_11"; 
  VariableNames[932] = "DPMS_DS_SURFACE_SPECIES_12"; 
  VariableNames[933] = "DPMS_DS_SURFACE_SPECIES_13"; 
  VariableNames[934] = "DPMS_DS_SURFACE_SPECIES_14"; 
  VariableNames[935] = "DPMS_DS_SURFACE_SPECIES_15"; 
  VariableNames[936] = "DPMS_DS_SURFACE_SPECIES_16"; 
  VariableNames[937] = "DPMS_DS_SURFACE_SPECIES_17"; 
  VariableNames[938] = "DPMS_DS_SURFACE_SPECIES_18"; 
  VariableNames[939] = "DPMS_DS_SURFACE_SPECIES_19"; 
  VariableNames[940] = "DPMS_DS_SURFACE_SPECIES_20"; 
  VariableNames[941] = "DPMS_DS_SURFACE_SPECIES_21"; 
  VariableNames[942] = "DPMS_DS_SURFACE_SPECIES_22"; 
  VariableNames[943] = "DPMS_DS_SURFACE_SPECIES_23"; 
  VariableNames[944] = "DPMS_DS_SURFACE_SPECIES_24"; 
  VariableNames[945] = "DPMS_DS_SURFACE_SPECIES_25"; 
  VariableNames[946] = "DPMS_DS_SURFACE_SPECIES_26"; 
  VariableNames[947] = "DPMS_DS_SURFACE_SPECIES_27"; 
  VariableNames[948] = "DPMS_DS_SURFACE_SPECIES_28"; 
  VariableNames[949] = "DPMS_DS_SURFACE_SPECIES_29"; 
  VariableNames[950] = "DPMS_DS_SURFACE_SPECIES_30"; 
  VariableNames[951] = "DPMS_DS_SURFACE_SPECIES_31"; 
  VariableNames[952] = "DPMS_DS_SURFACE_SPECIES_32"; 
  VariableNames[953] = "DPMS_DS_SURFACE_SPECIES_33"; 
  VariableNames[954] = "DPMS_DS_SURFACE_SPECIES_34"; 
  VariableNames[955] = "DPMS_DS_SURFACE_SPECIES_35"; 
  VariableNames[956] = "DPMS_DS_SURFACE_SPECIES_36"; 
  VariableNames[957] = "DPMS_DS_SURFACE_SPECIES_37"; 
  VariableNames[958] = "DPMS_DS_SURFACE_SPECIES_38"; 
  VariableNames[959] = "DPMS_DS_SURFACE_SPECIES_39"; 
  VariableNames[960] = "DPMS_DS_SURFACE_SPECIES_40"; 
  VariableNames[961] = "DPMS_DS_SURFACE_SPECIES_41"; 
  VariableNames[962] = "DPMS_DS_SURFACE_SPECIES_42"; 
  VariableNames[963] = "DPMS_DS_SURFACE_SPECIES_43"; 
  VariableNames[964] = "DPMS_DS_SURFACE_SPECIES_44"; 
  VariableNames[965] = "DPMS_DS_SURFACE_SPECIES_45"; 
  VariableNames[966] = "DPMS_DS_SURFACE_SPECIES_46"; 
  VariableNames[967] = "DPMS_DS_SURFACE_SPECIES_47"; 
  VariableNames[968] = "DPMS_DS_SURFACE_SPECIES_48"; 
  VariableNames[969] = "DPMS_DS_SURFACE_SPECIES_49";

  VariableNames[970] = "UDM_I";
 

  VariableNames[1000] = "Y_MEAN_00"; 
  VariableNames[1001] = "Y_MEAN_01"; 
  VariableNames[1002] = "Y_MEAN_02"; 
  VariableNames[1003] = "Y_MEAN_03"; 
  VariableNames[1004] = "Y_MEAN_04"; 
  VariableNames[1005] = "Y_MEAN_05"; 
  VariableNames[1006] = "Y_MEAN_06"; 
  VariableNames[1007] = "Y_MEAN_07"; 
  VariableNames[1008] = "Y_MEAN_08"; 
  VariableNames[1009] = "Y_MEAN_09"; 
  VariableNames[1010] = "Y_MEAN_10"; 
  VariableNames[1011] = "Y_MEAN_11"; 
  VariableNames[1012] = "Y_MEAN_12"; 
  VariableNames[1013] = "Y_MEAN_13"; 
  VariableNames[1014] = "Y_MEAN_14"; 
  VariableNames[1015] = "Y_MEAN_15"; 
  VariableNames[1016] = "Y_MEAN_16"; 
  VariableNames[1017] = "Y_MEAN_17"; 
  VariableNames[1018] = "Y_MEAN_18"; 
  VariableNames[1019] = "Y_MEAN_19"; 
  VariableNames[1020] = "Y_MEAN_20"; 
  VariableNames[1021] = "Y_MEAN_21"; 
  VariableNames[1022] = "Y_MEAN_22"; 
  VariableNames[1023] = "Y_MEAN_23"; 
  VariableNames[1024] = "Y_MEAN_24"; 
  VariableNames[1025] = "Y_MEAN_25"; 
  VariableNames[1026] = "Y_MEAN_26"; 
  VariableNames[1027] = "Y_MEAN_27"; 
  VariableNames[1028] = "Y_MEAN_28"; 
  VariableNames[1029] = "Y_MEAN_29"; 
  VariableNames[1030] = "Y_MEAN_30"; 
  VariableNames[1031] = "Y_MEAN_31"; 
  VariableNames[1032] = "Y_MEAN_32"; 
  VariableNames[1033] = "Y_MEAN_33"; 
  VariableNames[1034] = "Y_MEAN_34"; 
  VariableNames[1035] = "Y_MEAN_35"; 
  VariableNames[1036] = "Y_MEAN_36"; 
  VariableNames[1037] = "Y_MEAN_37"; 
  VariableNames[1038] = "Y_MEAN_38"; 
  VariableNames[1039] = "Y_MEAN_39"; 
  VariableNames[1040] = "Y_MEAN_40"; 
  VariableNames[1041] = "Y_MEAN_41"; 
  VariableNames[1042] = "Y_MEAN_42"; 
  VariableNames[1043] = "Y_MEAN_43"; 
  VariableNames[1044] = "Y_MEAN_44"; 
  VariableNames[1045] = "Y_MEAN_45"; 
  VariableNames[1046] = "Y_MEAN_46"; 
  VariableNames[1047] = "Y_MEAN_47"; 
  VariableNames[1048] = "Y_MEAN_48"; 
  VariableNames[1049] = "Y_MEAN_49"; 

  VariableNames[1050] = "Y_RMS_00"; 
  VariableNames[1051] = "Y_RMS_01"; 
  VariableNames[1052] = "Y_RMS_02"; 
  VariableNames[1053] = "Y_RMS_03"; 
  VariableNames[1054] = "Y_RMS_04"; 
  VariableNames[1055] = "Y_RMS_05"; 
  VariableNames[1056] = "Y_RMS_06"; 
  VariableNames[1057] = "Y_RMS_07"; 
  VariableNames[1058] = "Y_RMS_08"; 
  VariableNames[1059] = "Y_RMS_09"; 
  VariableNames[1060] = "Y_RMS_10"; 
  VariableNames[1061] = "Y_RMS_11"; 
  VariableNames[1062] = "Y_RMS_12"; 
  VariableNames[1063] = "Y_RMS_13"; 
  VariableNames[1064] = "Y_RMS_14"; 
  VariableNames[1065] = "Y_RMS_15"; 
  VariableNames[1066] = "Y_RMS_16"; 
  VariableNames[1067] = "Y_RMS_17"; 
  VariableNames[1068] = "Y_RMS_18"; 
  VariableNames[1069] = "Y_RMS_19"; 
  VariableNames[1070] = "Y_RMS_20"; 
  VariableNames[1071] = "Y_RMS_21"; 
  VariableNames[1072] = "Y_RMS_22"; 
  VariableNames[1073] = "Y_RMS_23"; 
  VariableNames[1074] = "Y_RMS_24"; 
  VariableNames[1075] = "Y_RMS_25"; 
  VariableNames[1076] = "Y_RMS_26"; 
  VariableNames[1077] = "Y_RMS_27"; 
  VariableNames[1078] = "Y_RMS_28"; 
  VariableNames[1079] = "Y_RMS_29"; 
  VariableNames[1080] = "Y_RMS_30"; 
  VariableNames[1081] = "Y_RMS_31"; 
  VariableNames[1082] = "Y_RMS_32"; 
  VariableNames[1083] = "Y_RMS_33"; 
  VariableNames[1084] = "Y_RMS_34"; 
  VariableNames[1085] = "Y_RMS_35"; 
  VariableNames[1086] = "Y_RMS_36"; 
  VariableNames[1087] = "Y_RMS_37"; 
  VariableNames[1088] = "Y_RMS_38"; 
  VariableNames[1089] = "Y_RMS_39"; 
  VariableNames[1090] = "Y_RMS_40"; 
  VariableNames[1091] = "Y_RMS_41"; 
  VariableNames[1092] = "Y_RMS_42"; 
  VariableNames[1093] = "Y_RMS_43"; 
  VariableNames[1094] = "Y_RMS_44"; 
  VariableNames[1095] = "Y_RMS_45"; 
  VariableNames[1096] = "Y_RMS_46"; 
  VariableNames[1097] = "Y_RMS_47"; 
  VariableNames[1098] = "Y_RMS_48"; 
  VariableNames[1099] = "Y_RMS_49"; 


  VariableNames[1200] = "SITE_F_00"; 
  VariableNames[1201] = "SITE_F_01"; 
  VariableNames[1202] = "SITE_F_02"; 
  VariableNames[1203] = "SITE_F_03"; 
  VariableNames[1204] = "SITE_F_04"; 
  VariableNames[1205] = "SITE_F_05"; 
  VariableNames[1206] = "SITE_F_06"; 
  VariableNames[1207] = "SITE_F_07"; 
  VariableNames[1208] = "SITE_F_08"; 
  VariableNames[1209] = "SITE_F_09"; 
  VariableNames[1210] = "SITE_F_10"; 
  VariableNames[1211] = "SITE_F_11"; 
  VariableNames[1212] = "SITE_F_12"; 
  VariableNames[1213] = "SITE_F_13"; 
  VariableNames[1214] = "SITE_F_14"; 
  VariableNames[1215] = "SITE_F_15"; 
  VariableNames[1216] = "SITE_F_16"; 
  VariableNames[1217] = "SITE_F_17"; 
  VariableNames[1218] = "SITE_F_18"; 
  VariableNames[1219] = "SITE_F_19"; 
  VariableNames[1220] = "SITE_F_20"; 
  VariableNames[1221] = "SITE_F_21"; 
  VariableNames[1222] = "SITE_F_22"; 
  VariableNames[1223] = "SITE_F_23"; 
  VariableNames[1224] = "SITE_F_24"; 
  VariableNames[1225] = "SITE_F_25"; 
  VariableNames[1226] = "SITE_F_26"; 
  VariableNames[1227] = "SITE_F_27"; 
  VariableNames[1228] = "SITE_F_28"; 
  VariableNames[1229] = "SITE_F_29"; 
  VariableNames[1230] = "SITE_F_30"; 
  VariableNames[1231] = "SITE_F_31"; 
  VariableNames[1232] = "SITE_F_32"; 
  VariableNames[1233] = "SITE_F_33"; 
  VariableNames[1234] = "SITE_F_34"; 
  VariableNames[1235] = "SITE_F_35"; 
  VariableNames[1236] = "SITE_F_36"; 
  VariableNames[1237] = "SITE_F_37"; 
  VariableNames[1238] = "SITE_F_38"; 
  VariableNames[1239] = "SITE_F_39"; 
  VariableNames[1240] = "SITE_F_40"; 
  VariableNames[1241] = "SITE_F_41"; 
  VariableNames[1242] = "SITE_F_42"; 
  VariableNames[1243] = "SITE_F_43"; 
  VariableNames[1244] = "SITE_F_44"; 
  VariableNames[1245] = "SITE_F_45"; 
  VariableNames[1246] = "SITE_F_46"; 
  VariableNames[1247] = "SITE_F_47"; 
  VariableNames[1248] = "SITE_F_48"; 
  VariableNames[1249] = "SITE_F_49"; 

  VariableNames[1250] = "CREV_Y_00"; 
  VariableNames[1251] = "CREV_Y_01"; 
  VariableNames[1252] = "CREV_Y_02"; 
  VariableNames[1253] = "CREV_Y_03"; 
  VariableNames[1254] = "CREV_Y_04"; 
  VariableNames[1255] = "CREV_Y_05"; 
  VariableNames[1256] = "CREV_Y_06"; 
  VariableNames[1257] = "CREV_Y_07"; 
  VariableNames[1258] = "CREV_Y_08"; 
  VariableNames[1259] = "CREV_Y_09"; 
  VariableNames[1260] = "CREV_Y_10"; 
  VariableNames[1261] = "CREV_Y_11"; 
  VariableNames[1262] = "CREV_Y_12"; 
  VariableNames[1263] = "CREV_Y_13"; 
  VariableNames[1264] = "CREV_Y_14"; 
  VariableNames[1265] = "CREV_Y_15"; 
  VariableNames[1266] = "CREV_Y_16"; 
  VariableNames[1267] = "CREV_Y_17"; 
  VariableNames[1268] = "CREV_Y_18"; 
  VariableNames[1269] = "CREV_Y_19"; 
  VariableNames[1270] = "CREV_Y_20"; 
  VariableNames[1271] = "CREV_Y_21"; 
  VariableNames[1272] = "CREV_Y_22"; 
  VariableNames[1273] = "CREV_Y_23"; 
  VariableNames[1274] = "CREV_Y_24"; 
  VariableNames[1275] = "CREV_Y_25"; 
  VariableNames[1276] = "CREV_Y_26"; 
  VariableNames[1277] = "CREV_Y_27"; 
  VariableNames[1278] = "CREV_Y_28"; 
  VariableNames[1279] = "CREV_Y_29"; 
  VariableNames[1280] = "CREV_Y_30"; 
  VariableNames[1281] = "CREV_Y_31"; 
  VariableNames[1282] = "CREV_Y_32"; 
  VariableNames[1283] = "CREV_Y_33"; 
  VariableNames[1284] = "CREV_Y_34"; 
  VariableNames[1285] = "CREV_Y_35"; 
  VariableNames[1286] = "CREV_Y_36"; 
  VariableNames[1287] = "CREV_Y_37"; 
  VariableNames[1288] = "CREV_Y_38"; 
  VariableNames[1289] = "CREV_Y_39"; 
  VariableNames[1290] = "CREV_Y_40"; 
  VariableNames[1291] = "CREV_Y_41"; 
  VariableNames[1292] = "CREV_Y_42"; 
  VariableNames[1293] = "CREV_Y_43"; 
  VariableNames[1294] = "CREV_Y_44"; 
  VariableNames[1295] = "CREV_Y_45"; 
  VariableNames[1296] = "CREV_Y_46"; 
  VariableNames[1297] = "CREV_Y_47"; 
  VariableNames[1298] = "CREV_Y_48"; 
  VariableNames[1299] = "CREV_Y_49"; 

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

  VariableNames[1350] = "SURF_F_00"; 
  VariableNames[1351] = "SURF_F_01"; 
  VariableNames[1352] = "SURF_F_02"; 
  VariableNames[1353] = "SURF_F_03"; 
  VariableNames[1354] = "SURF_F_04"; 
  VariableNames[1355] = "SURF_F_05"; 
  VariableNames[1356] = "SURF_F_06"; 
  VariableNames[1357] = "SURF_F_07"; 
  VariableNames[1358] = "SURF_F_08"; 
  VariableNames[1359] = "SURF_F_09"; 
  VariableNames[1360] = "SURF_F_10"; 
  VariableNames[1361] = "SURF_F_11"; 
  VariableNames[1362] = "SURF_F_12"; 
  VariableNames[1363] = "SURF_F_13"; 
  VariableNames[1364] = "SURF_F_14"; 
  VariableNames[1365] = "SURF_F_15"; 
  VariableNames[1366] = "SURF_F_16"; 
  VariableNames[1367] = "SURF_F_17"; 
  VariableNames[1368] = "SURF_F_18"; 
  VariableNames[1369] = "SURF_F_19"; 
  VariableNames[1370] = "SURF_F_20"; 
  VariableNames[1371] = "SURF_F_21"; 
  VariableNames[1372] = "SURF_F_22"; 
  VariableNames[1373] = "SURF_F_23"; 
  VariableNames[1374] = "SURF_F_24"; 
  VariableNames[1375] = "SURF_F_25"; 
  VariableNames[1376] = "SURF_F_26"; 
  VariableNames[1377] = "SURF_F_27"; 
  VariableNames[1378] = "SURF_F_28"; 
  VariableNames[1379] = "SURF_F_29"; 
  VariableNames[1380] = "SURF_F_30"; 
  VariableNames[1381] = "SURF_F_31"; 
  VariableNames[1382] = "SURF_F_32"; 
  VariableNames[1383] = "SURF_F_33"; 
  VariableNames[1384] = "SURF_F_34"; 
  VariableNames[1385] = "SURF_F_35"; 
  VariableNames[1386] = "SURF_F_36"; 
  VariableNames[1387] = "SURF_F_37"; 
  VariableNames[1388] = "SURF_F_38"; 
  VariableNames[1389] = "SURF_F_39"; 
  VariableNames[1390] = "SURF_F_40"; 
  VariableNames[1391] = "SURF_F_41"; 
  VariableNames[1392] = "SURF_F_42"; 
  VariableNames[1393] = "SURF_F_43"; 
  VariableNames[1394] = "SURF_F_44"; 
  VariableNames[1395] = "SURF_F_45"; 
  VariableNames[1396] = "SURF_F_46"; 
  VariableNames[1397] = "SURF_F_47"; 
  VariableNames[1398] = "SURF_F_48"; 
  VariableNames[1399] = "SURF_F_49"; 








  VariableNames[7700] = "PB_DISC_00"; 
  VariableNames[7701] = "PB_DISC_01"; 
  VariableNames[7702] = "PB_DISC_02"; 
  VariableNames[7703] = "PB_DISC_03"; 
  VariableNames[7704] = "PB_DISC_04"; 
  VariableNames[7705] = "PB_DISC_05"; 
  VariableNames[7706] = "PB_DISC_06"; 
  VariableNames[7707] = "PB_DISC_07"; 
  VariableNames[7708] = "PB_DISC_08"; 
  VariableNames[7709] = "PB_DISC_09"; 
  VariableNames[7710] = "PB_DISC_10"; 
  VariableNames[7711] = "PB_DISC_11"; 
  VariableNames[7712] = "PB_DISC_12"; 
  VariableNames[7713] = "PB_DISC_13"; 
  VariableNames[7714] = "PB_DISC_14"; 
  VariableNames[7715] = "PB_DISC_15"; 
  VariableNames[7716] = "PB_DISC_16"; 
  VariableNames[7717] = "PB_DISC_17"; 
  VariableNames[7718] = "PB_DISC_18"; 
  VariableNames[7719] = "PB_DISC_19"; 
  VariableNames[7720] = "PB_DISC_20"; 
  VariableNames[7721] = "PB_DISC_21"; 
  VariableNames[7722] = "PB_DISC_22"; 
  VariableNames[7723] = "PB_DISC_23"; 
  VariableNames[7724] = "PB_DISC_24"; 
  VariableNames[7725] = "PB_DISC_25"; 
  VariableNames[7726] = "PB_DISC_26"; 
  VariableNames[7727] = "PB_DISC_27"; 
  VariableNames[7728] = "PB_DISC_28"; 
  VariableNames[7729] = "PB_DISC_29"; 
  VariableNames[7730] = "PB_DISC_30"; 
  VariableNames[7731] = "PB_DISC_31"; 
  VariableNames[7732] = "PB_DISC_32"; 
  VariableNames[7733] = "PB_DISC_33"; 
  VariableNames[7734] = "PB_DISC_34"; 
  VariableNames[7735] = "PB_DISC_35"; 
  VariableNames[7736] = "PB_DISC_36"; 
  VariableNames[7737] = "PB_DISC_37"; 
  VariableNames[7738] = "PB_DISC_38"; 
  VariableNames[7739] = "PB_DISC_39"; 
  VariableNames[7740] = "PB_DISC_40"; 
  VariableNames[7741] = "PB_DISC_41"; 
  VariableNames[7742] = "PB_DISC_42"; 
  VariableNames[7743] = "PB_DISC_43"; 
  VariableNames[7744] = "PB_DISC_44"; 
  VariableNames[7745] = "PB_DISC_45"; 
  VariableNames[7746] = "PB_DISC_46"; 
  VariableNames[7747] = "PB_DISC_47"; 
  VariableNames[7748] = "PB_DISC_48"; 
  VariableNames[7749] = "PB_DISC_49"; 

  VariableNames[7750] = "PB_DISC_M1_00"; 
  VariableNames[7751] = "PB_DISC_M1_01"; 
  VariableNames[7752] = "PB_DISC_M1_02"; 
  VariableNames[7753] = "PB_DISC_M1_03"; 
  VariableNames[7754] = "PB_DISC_M1_04"; 
  VariableNames[7755] = "PB_DISC_M1_05"; 
  VariableNames[7756] = "PB_DISC_M1_06"; 
  VariableNames[7757] = "PB_DISC_M1_07"; 
  VariableNames[7758] = "PB_DISC_M1_08"; 
  VariableNames[7759] = "PB_DISC_M1_09"; 
  VariableNames[7760] = "PB_DISC_M1_10"; 
  VariableNames[7761] = "PB_DISC_M1_11"; 
  VariableNames[7762] = "PB_DISC_M1_12"; 
  VariableNames[7763] = "PB_DISC_M1_13"; 
  VariableNames[7764] = "PB_DISC_M1_14"; 
  VariableNames[7765] = "PB_DISC_M1_15"; 
  VariableNames[7766] = "PB_DISC_M1_16"; 
  VariableNames[7767] = "PB_DISC_M1_17"; 
  VariableNames[7768] = "PB_DISC_M1_18"; 
  VariableNames[7769] = "PB_DISC_M1_19"; 
  VariableNames[7770] = "PB_DISC_M1_20"; 
  VariableNames[7771] = "PB_DISC_M1_21"; 
  VariableNames[7772] = "PB_DISC_M1_22"; 
  VariableNames[7773] = "PB_DISC_M1_23"; 
  VariableNames[7774] = "PB_DISC_M1_24"; 
  VariableNames[7775] = "PB_DISC_M1_25"; 
  VariableNames[7776] = "PB_DISC_M1_26"; 
  VariableNames[7777] = "PB_DISC_M1_27"; 
  VariableNames[7778] = "PB_DISC_M1_28"; 
  VariableNames[7779] = "PB_DISC_M1_29"; 
  VariableNames[7780] = "PB_DISC_M1_30"; 
  VariableNames[7781] = "PB_DISC_M1_31"; 
  VariableNames[7782] = "PB_DISC_M1_32"; 
  VariableNames[7783] = "PB_DISC_M1_33"; 
  VariableNames[7784] = "PB_DISC_M1_34"; 
  VariableNames[7785] = "PB_DISC_M1_35"; 
  VariableNames[7786] = "PB_DISC_M1_36"; 
  VariableNames[7787] = "PB_DISC_M1_37"; 
  VariableNames[7788] = "PB_DISC_M1_38"; 
  VariableNames[7789] = "PB_DISC_M1_39"; 
  VariableNames[7790] = "PB_DISC_M1_40"; 
  VariableNames[7791] = "PB_DISC_M1_41"; 
  VariableNames[7792] = "PB_DISC_M1_42"; 
  VariableNames[7793] = "PB_DISC_M1_43"; 
  VariableNames[7794] = "PB_DISC_M1_44"; 
  VariableNames[7795] = "PB_DISC_M1_45"; 
  VariableNames[7796] = "PB_DISC_M1_46"; 
  VariableNames[7797] = "PB_DISC_M1_47"; 
  VariableNames[7798] = "PB_DISC_M1_48"; 
  VariableNames[7799] = "PB_DISC_M1_49"; 

  VariableNames[7800] = "PB_DISC_M2_00"; 
  VariableNames[7801] = "PB_DISC_M2_01"; 
  VariableNames[7802] = "PB_DISC_M2_02"; 
  VariableNames[7803] = "PB_DISC_M2_03"; 
  VariableNames[7804] = "PB_DISC_M2_04"; 
  VariableNames[7805] = "PB_DISC_M2_05"; 
  VariableNames[7806] = "PB_DISC_M2_06"; 
  VariableNames[7807] = "PB_DISC_M2_07"; 
  VariableNames[7808] = "PB_DISC_M2_08"; 
  VariableNames[7809] = "PB_DISC_M2_09"; 
  VariableNames[7810] = "PB_DISC_M2_10"; 
  VariableNames[7811] = "PB_DISC_M2_11"; 
  VariableNames[7812] = "PB_DISC_M2_12"; 
  VariableNames[7813] = "PB_DISC_M2_13"; 
  VariableNames[7814] = "PB_DISC_M2_14"; 
  VariableNames[7815] = "PB_DISC_M2_15"; 
  VariableNames[7816] = "PB_DISC_M2_16"; 
  VariableNames[7817] = "PB_DISC_M2_17"; 
  VariableNames[7818] = "PB_DISC_M2_18"; 
  VariableNames[7819] = "PB_DISC_M2_19"; 
  VariableNames[7820] = "PB_DISC_M2_20"; 
  VariableNames[7821] = "PB_DISC_M2_21"; 
  VariableNames[7822] = "PB_DISC_M2_22"; 
  VariableNames[7823] = "PB_DISC_M2_23"; 
  VariableNames[7824] = "PB_DISC_M2_24"; 
  VariableNames[7825] = "PB_DISC_M2_25"; 
  VariableNames[7826] = "PB_DISC_M2_26"; 
  VariableNames[7827] = "PB_DISC_M2_27"; 
  VariableNames[7828] = "PB_DISC_M2_28"; 
  VariableNames[7829] = "PB_DISC_M2_29"; 
  VariableNames[7830] = "PB_DISC_M2_30"; 
  VariableNames[7831] = "PB_DISC_M2_31"; 
  VariableNames[7832] = "PB_DISC_M2_32"; 
  VariableNames[7833] = "PB_DISC_M2_33"; 
  VariableNames[7834] = "PB_DISC_M2_34"; 
  VariableNames[7835] = "PB_DISC_M2_35"; 
  VariableNames[7836] = "PB_DISC_M2_36"; 
  VariableNames[7837] = "PB_DISC_M2_37"; 
  VariableNames[7838] = "PB_DISC_M2_38"; 
  VariableNames[7839] = "PB_DISC_M2_39"; 
  VariableNames[7840] = "PB_DISC_M2_40"; 
  VariableNames[7841] = "PB_DISC_M2_41"; 
  VariableNames[7842] = "PB_DISC_M2_42"; 
  VariableNames[7843] = "PB_DISC_M2_43"; 
  VariableNames[7844] = "PB_DISC_M2_44"; 
  VariableNames[7845] = "PB_DISC_M2_45"; 
  VariableNames[7846] = "PB_DISC_M2_46"; 
  VariableNames[7847] = "PB_DISC_M2_47"; 
  VariableNames[7848] = "PB_DISC_M2_48"; 
  VariableNames[7849] = "PB_DISC_M2_49"; 

  VariableNames[7850] = "PB_QMOM_00"; 
  VariableNames[7851] = "PB_QMOM_01"; 
  VariableNames[7852] = "PB_QMOM_02"; 
  VariableNames[7853] = "PB_QMOM_03"; 
  VariableNames[7854] = "PB_QMOM_04"; 
  VariableNames[7855] = "PB_QMOM_05"; 
  VariableNames[7856] = "PB_QMOM_06"; 
  VariableNames[7857] = "PB_QMOM_07"; 
  VariableNames[7858] = "PB_QMOM_08"; 
  VariableNames[7859] = "PB_QMOM_09"; 
  VariableNames[7860] = "PB_QMOM_10"; 
  VariableNames[7861] = "PB_QMOM_11"; 
  VariableNames[7862] = "PB_QMOM_12"; 
  VariableNames[7863] = "PB_QMOM_13"; 
  VariableNames[7864] = "PB_QMOM_14"; 
  VariableNames[7865] = "PB_QMOM_15"; 
  VariableNames[7866] = "PB_QMOM_16"; 
  VariableNames[7867] = "PB_QMOM_17"; 
  VariableNames[7868] = "PB_QMOM_18"; 
  VariableNames[7869] = "PB_QMOM_19"; 
  VariableNames[7870] = "PB_QMOM_20"; 
  VariableNames[7871] = "PB_QMOM_21"; 
  VariableNames[7872] = "PB_QMOM_22"; 
  VariableNames[7873] = "PB_QMOM_23"; 
  VariableNames[7874] = "PB_QMOM_24"; 
  VariableNames[7875] = "PB_QMOM_25"; 
  VariableNames[7876] = "PB_QMOM_26"; 
  VariableNames[7877] = "PB_QMOM_27"; 
  VariableNames[7878] = "PB_QMOM_28"; 
  VariableNames[7879] = "PB_QMOM_29"; 
  VariableNames[7880] = "PB_QMOM_30"; 
  VariableNames[7881] = "PB_QMOM_31"; 
  VariableNames[7882] = "PB_QMOM_32"; 
  VariableNames[7883] = "PB_QMOM_33"; 
  VariableNames[7884] = "PB_QMOM_34"; 
  VariableNames[7885] = "PB_QMOM_35"; 
  VariableNames[7886] = "PB_QMOM_36"; 
  VariableNames[7887] = "PB_QMOM_37"; 
  VariableNames[7888] = "PB_QMOM_38"; 
  VariableNames[7889] = "PB_QMOM_39"; 
  VariableNames[7890] = "PB_QMOM_40"; 
  VariableNames[7891] = "PB_QMOM_41"; 
  VariableNames[7892] = "PB_QMOM_42"; 
  VariableNames[7893] = "PB_QMOM_43"; 
  VariableNames[7894] = "PB_QMOM_44"; 
  VariableNames[7895] = "PB_QMOM_45"; 
  VariableNames[7896] = "PB_QMOM_46"; 
  VariableNames[7897] = "PB_QMOM_47"; 
  VariableNames[7898] = "PB_QMOM_48"; 
  VariableNames[7899] = "PB_QMOM_49"; 

  VariableNames[7900] = "PB_QMOM_M1_00"; 
  VariableNames[7901] = "PB_QMOM_M1_01"; 
  VariableNames[7902] = "PB_QMOM_M1_02"; 
  VariableNames[7903] = "PB_QMOM_M1_03"; 
  VariableNames[7904] = "PB_QMOM_M1_04"; 
  VariableNames[7905] = "PB_QMOM_M1_05"; 
  VariableNames[7906] = "PB_QMOM_M1_06"; 
  VariableNames[7907] = "PB_QMOM_M1_07"; 
  VariableNames[7908] = "PB_QMOM_M1_08"; 
  VariableNames[7909] = "PB_QMOM_M1_09"; 
  VariableNames[7910] = "PB_QMOM_M1_10"; 
  VariableNames[7911] = "PB_QMOM_M1_11"; 
  VariableNames[7912] = "PB_QMOM_M1_12"; 
  VariableNames[7913] = "PB_QMOM_M1_13"; 
  VariableNames[7914] = "PB_QMOM_M1_14"; 
  VariableNames[7915] = "PB_QMOM_M1_15"; 
  VariableNames[7916] = "PB_QMOM_M1_16"; 
  VariableNames[7917] = "PB_QMOM_M1_17"; 
  VariableNames[7918] = "PB_QMOM_M1_18"; 
  VariableNames[7919] = "PB_QMOM_M1_19"; 
  VariableNames[7920] = "PB_QMOM_M1_20"; 
  VariableNames[7921] = "PB_QMOM_M1_21"; 
  VariableNames[7922] = "PB_QMOM_M1_22"; 
  VariableNames[7923] = "PB_QMOM_M1_23"; 
  VariableNames[7924] = "PB_QMOM_M1_24"; 
  VariableNames[7925] = "PB_QMOM_M1_25"; 
  VariableNames[7926] = "PB_QMOM_M1_26"; 
  VariableNames[7927] = "PB_QMOM_M1_27"; 
  VariableNames[7928] = "PB_QMOM_M1_28"; 
  VariableNames[7929] = "PB_QMOM_M1_29"; 
  VariableNames[7930] = "PB_QMOM_M1_30"; 
  VariableNames[7931] = "PB_QMOM_M1_31"; 
  VariableNames[7932] = "PB_QMOM_M1_32"; 
  VariableNames[7933] = "PB_QMOM_M1_33"; 
  VariableNames[7934] = "PB_QMOM_M1_34"; 
  VariableNames[7935] = "PB_QMOM_M1_35"; 
  VariableNames[7936] = "PB_QMOM_M1_36"; 
  VariableNames[7937] = "PB_QMOM_M1_37"; 
  VariableNames[7938] = "PB_QMOM_M1_38"; 
  VariableNames[7939] = "PB_QMOM_M1_39"; 
  VariableNames[7940] = "PB_QMOM_M1_40"; 
  VariableNames[7941] = "PB_QMOM_M1_41"; 
  VariableNames[7942] = "PB_QMOM_M1_42"; 
  VariableNames[7943] = "PB_QMOM_M1_43"; 
  VariableNames[7944] = "PB_QMOM_M1_44"; 
  VariableNames[7945] = "PB_QMOM_M1_45"; 
  VariableNames[7946] = "PB_QMOM_M1_46"; 
  VariableNames[7947] = "PB_QMOM_M1_47"; 
  VariableNames[7948] = "PB_QMOM_M1_48"; 
  VariableNames[7949] = "PB_QMOM_M1_49"; 


  VariableNames[7950] = "PB_QMOM_M2_00"; 
  VariableNames[7951] = "PB_QMOM_M2_01"; 
  VariableNames[7952] = "PB_QMOM_M2_02"; 
  VariableNames[7953] = "PB_QMOM_M2_03"; 
  VariableNames[7954] = "PB_QMOM_M2_04"; 
  VariableNames[7955] = "PB_QMOM_M2_05"; 
  VariableNames[7956] = "PB_QMOM_M2_06"; 
  VariableNames[7957] = "PB_QMOM_M2_07"; 
  VariableNames[7958] = "PB_QMOM_M2_08"; 
  VariableNames[7959] = "PB_QMOM_M2_09"; 
  VariableNames[7960] = "PB_QMOM_M2_10"; 
  VariableNames[7961] = "PB_QMOM_M2_11"; 
  VariableNames[7962] = "PB_QMOM_M2_12"; 
  VariableNames[7963] = "PB_QMOM_M2_13"; 
  VariableNames[7964] = "PB_QMOM_M2_14"; 
  VariableNames[7965] = "PB_QMOM_M2_15"; 
  VariableNames[7966] = "PB_QMOM_M2_16"; 
  VariableNames[7967] = "PB_QMOM_M2_17"; 
  VariableNames[7968] = "PB_QMOM_M2_18"; 
  VariableNames[7969] = "PB_QMOM_M2_19"; 
  VariableNames[7970] = "PB_QMOM_M2_20"; 
  VariableNames[7971] = "PB_QMOM_M2_21"; 
  VariableNames[7972] = "PB_QMOM_M2_22"; 
  VariableNames[7973] = "PB_QMOM_M2_23"; 
  VariableNames[7974] = "PB_QMOM_M2_24"; 
  VariableNames[7975] = "PB_QMOM_M2_25"; 
  VariableNames[7976] = "PB_QMOM_M2_26"; 
  VariableNames[7977] = "PB_QMOM_M2_27"; 
  VariableNames[7978] = "PB_QMOM_M2_28"; 
  VariableNames[7979] = "PB_QMOM_M2_29"; 
  VariableNames[7980] = "PB_QMOM_M2_30"; 
  VariableNames[7981] = "PB_QMOM_M2_31"; 
  VariableNames[7982] = "PB_QMOM_M2_32"; 
  VariableNames[7983] = "PB_QMOM_M2_33"; 
  VariableNames[7984] = "PB_QMOM_M2_34"; 
  VariableNames[7985] = "PB_QMOM_M2_35"; 
  VariableNames[7986] = "PB_QMOM_M2_36"; 
  VariableNames[7987] = "PB_QMOM_M2_37"; 
  VariableNames[7988] = "PB_QMOM_M2_38"; 
  VariableNames[7989] = "PB_QMOM_M2_39"; 
  VariableNames[7990] = "PB_QMOM_M2_40"; 
  VariableNames[7991] = "PB_QMOM_M2_41"; 
  VariableNames[7992] = "PB_QMOM_M2_42"; 
  VariableNames[7993] = "PB_QMOM_M2_43"; 
  VariableNames[7994] = "PB_QMOM_M2_44"; 
  VariableNames[7995] = "PB_QMOM_M2_45"; 
  VariableNames[7996] = "PB_QMOM_M2_46"; 
  VariableNames[7997] = "PB_QMOM_M2_47"; 
  VariableNames[7998] = "PB_QMOM_M2_48"; 
  VariableNames[7999] = "PB_QMOM_M2_49"; 

  VariableNames[8000] = "PB_SMM_00"; 
  VariableNames[8001] = "PB_SMM_01"; 
  VariableNames[8002] = "PB_SMM_02"; 
  VariableNames[8003] = "PB_SMM_03"; 
  VariableNames[8004] = "PB_SMM_04"; 
  VariableNames[8005] = "PB_SMM_05"; 
  VariableNames[8006] = "PB_SMM_06"; 
  VariableNames[8007] = "PB_SMM_07"; 
  VariableNames[8008] = "PB_SMM_08"; 
  VariableNames[8009] = "PB_SMM_09"; 
  VariableNames[8010] = "PB_SMM_10"; 
  VariableNames[8011] = "PB_SMM_11"; 
  VariableNames[8012] = "PB_SMM_12"; 
  VariableNames[8013] = "PB_SMM_13"; 
  VariableNames[8014] = "PB_SMM_14"; 
  VariableNames[8015] = "PB_SMM_15"; 
  VariableNames[8016] = "PB_SMM_16"; 
  VariableNames[8017] = "PB_SMM_17"; 
  VariableNames[8018] = "PB_SMM_18"; 
  VariableNames[8019] = "PB_SMM_19"; 
  VariableNames[8020] = "PB_SMM_20"; 
  VariableNames[8021] = "PB_SMM_21"; 
  VariableNames[8022] = "PB_SMM_22"; 
  VariableNames[8023] = "PB_SMM_23"; 
  VariableNames[8024] = "PB_SMM_24"; 
  VariableNames[8025] = "PB_SMM_25"; 
  VariableNames[8026] = "PB_SMM_26"; 
  VariableNames[8027] = "PB_SMM_27"; 
  VariableNames[8028] = "PB_SMM_28"; 
  VariableNames[8029] = "PB_SMM_29"; 
  VariableNames[8030] = "PB_SMM_30"; 
  VariableNames[8031] = "PB_SMM_31"; 
  VariableNames[8032] = "PB_SMM_32"; 
  VariableNames[8033] = "PB_SMM_33"; 
  VariableNames[8034] = "PB_SMM_34"; 
  VariableNames[8035] = "PB_SMM_35"; 
  VariableNames[8036] = "PB_SMM_36"; 
  VariableNames[8037] = "PB_SMM_37"; 
  VariableNames[8038] = "PB_SMM_38"; 
  VariableNames[8039] = "PB_SMM_39"; 
  VariableNames[8040] = "PB_SMM_40"; 
  VariableNames[8041] = "PB_SMM_41"; 
  VariableNames[8042] = "PB_SMM_42"; 
  VariableNames[8043] = "PB_SMM_43"; 
  VariableNames[8044] = "PB_SMM_44"; 
  VariableNames[8045] = "PB_SMM_45"; 
  VariableNames[8046] = "PB_SMM_46"; 
  VariableNames[8047] = "PB_SMM_47"; 
  VariableNames[8048] = "PB_SMM_48"; 
  VariableNames[8049] = "PB_SMM_49"; 


  VariableNames[8050] = "PB_SMM_M1_00"; 
  VariableNames[8051] = "PB_SMM_M1_01"; 
  VariableNames[8052] = "PB_SMM_M1_02"; 
  VariableNames[8053] = "PB_SMM_M1_03"; 
  VariableNames[8054] = "PB_SMM_M1_04"; 
  VariableNames[8055] = "PB_SMM_M1_05"; 
  VariableNames[8056] = "PB_SMM_M1_06"; 
  VariableNames[8057] = "PB_SMM_M1_07"; 
  VariableNames[8058] = "PB_SMM_M1_08"; 
  VariableNames[8059] = "PB_SMM_M1_09"; 
  VariableNames[8060] = "PB_SMM_M1_10"; 
  VariableNames[8061] = "PB_SMM_M1_11"; 
  VariableNames[8062] = "PB_SMM_M1_12"; 
  VariableNames[8063] = "PB_SMM_M1_13"; 
  VariableNames[8064] = "PB_SMM_M1_14"; 
  VariableNames[8065] = "PB_SMM_M1_15"; 
  VariableNames[8066] = "PB_SMM_M1_16"; 
  VariableNames[8067] = "PB_SMM_M1_17"; 
  VariableNames[8068] = "PB_SMM_M1_18"; 
  VariableNames[8069] = "PB_SMM_M1_19"; 
  VariableNames[8070] = "PB_SMM_M1_20"; 
  VariableNames[8071] = "PB_SMM_M1_21"; 
  VariableNames[8072] = "PB_SMM_M1_22"; 
  VariableNames[8073] = "PB_SMM_M1_23"; 
  VariableNames[8074] = "PB_SMM_M1_24"; 
  VariableNames[8075] = "PB_SMM_M1_25"; 
  VariableNames[8076] = "PB_SMM_M1_26"; 
  VariableNames[8077] = "PB_SMM_M1_27"; 
  VariableNames[8078] = "PB_SMM_M1_28"; 
  VariableNames[8079] = "PB_SMM_M1_29"; 
  VariableNames[8080] = "PB_SMM_M1_30"; 
  VariableNames[8081] = "PB_SMM_M1_31"; 
  VariableNames[8082] = "PB_SMM_M1_32"; 
  VariableNames[8083] = "PB_SMM_M1_33"; 
  VariableNames[8084] = "PB_SMM_M1_34"; 
  VariableNames[8085] = "PB_SMM_M1_35"; 
  VariableNames[8086] = "PB_SMM_M1_36"; 
  VariableNames[8087] = "PB_SMM_M1_37"; 
  VariableNames[8088] = "PB_SMM_M1_38"; 
  VariableNames[8089] = "PB_SMM_M1_39"; 
  VariableNames[8090] = "PB_SMM_M1_40"; 
  VariableNames[8091] = "PB_SMM_M1_41"; 
  VariableNames[8092] = "PB_SMM_M1_42"; 
  VariableNames[8093] = "PB_SMM_M1_43"; 
  VariableNames[8094] = "PB_SMM_M1_44"; 
  VariableNames[8095] = "PB_SMM_M1_45"; 
  VariableNames[8096] = "PB_SMM_M1_46"; 
  VariableNames[8097] = "PB_SMM_M1_47"; 
  VariableNames[8098] = "PB_SMM_M1_48"; 
  VariableNames[8099] = "PB_SMM_M1_49"; 

  VariableNames[8100] = "PB_SMM_M2_00"; 
  VariableNames[8101] = "PB_SMM_M2_01"; 
  VariableNames[8102] = "PB_SMM_M2_02"; 
  VariableNames[8103] = "PB_SMM_M2_03"; 
  VariableNames[8104] = "PB_SMM_M2_04"; 
  VariableNames[8105] = "PB_SMM_M2_05"; 
  VariableNames[8106] = "PB_SMM_M2_06"; 
  VariableNames[8107] = "PB_SMM_M2_07"; 
  VariableNames[8108] = "PB_SMM_M2_08"; 
  VariableNames[8109] = "PB_SMM_M2_09"; 
  VariableNames[8110] = "PB_SMM_M2_10"; 
  VariableNames[8111] = "PB_SMM_M2_11"; 
  VariableNames[8112] = "PB_SMM_M2_12"; 
  VariableNames[8113] = "PB_SMM_M2_13"; 
  VariableNames[8114] = "PB_SMM_M2_14"; 
  VariableNames[8115] = "PB_SMM_M2_15"; 
  VariableNames[8116] = "PB_SMM_M2_16"; 
  VariableNames[8117] = "PB_SMM_M2_17"; 
  VariableNames[8118] = "PB_SMM_M2_18"; 
  VariableNames[8119] = "PB_SMM_M2_19"; 
  VariableNames[8120] = "PB_SMM_M2_20"; 
  VariableNames[8121] = "PB_SMM_M2_21"; 
  VariableNames[8122] = "PB_SMM_M2_22"; 
  VariableNames[8123] = "PB_SMM_M2_23"; 
  VariableNames[8124] = "PB_SMM_M2_24"; 
  VariableNames[8125] = "PB_SMM_M2_25"; 
  VariableNames[8126] = "PB_SMM_M2_26"; 
  VariableNames[8127] = "PB_SMM_M2_27"; 
  VariableNames[8128] = "PB_SMM_M2_28"; 
  VariableNames[8129] = "PB_SMM_M2_29"; 
  VariableNames[8130] = "PB_SMM_M2_30"; 
  VariableNames[8131] = "PB_SMM_M2_31"; 
  VariableNames[8132] = "PB_SMM_M2_32"; 
  VariableNames[8133] = "PB_SMM_M2_33"; 
  VariableNames[8134] = "PB_SMM_M2_34"; 
  VariableNames[8135] = "PB_SMM_M2_35"; 
  VariableNames[8136] = "PB_SMM_M2_36"; 
  VariableNames[8137] = "PB_SMM_M2_37"; 
  VariableNames[8138] = "PB_SMM_M2_38"; 
  VariableNames[8139] = "PB_SMM_M2_39"; 
  VariableNames[8140] = "PB_SMM_M2_40"; 
  VariableNames[8141] = "PB_SMM_M2_41"; 
  VariableNames[8142] = "PB_SMM_M2_42"; 
  VariableNames[8143] = "PB_SMM_M2_43"; 
  VariableNames[8144] = "PB_SMM_M2_44"; 
  VariableNames[8145] = "PB_SMM_M2_45"; 
  VariableNames[8146] = "PB_SMM_M2_46"; 
  VariableNames[8147] = "PB_SMM_M2_47"; 
  VariableNames[8148] = "PB_SMM_M2_48"; 
  VariableNames[8149] = "PB_SMM_M2_49"; 

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::ParseCaseFile()
{
  FluentCaseFile.clear();
  FluentCaseFile.seekg (0, ios::beg);

  int cnt = 0;
  while (GetCaseChunk())
    {

    int index = GetCaseIndex();
cout<<"parseCase: "<<index<<endl;
    switch (index)
      {
      case 0:
        break;
      case 1:
        break;
      case 2:
        GridDimension = GetDimension();
        break;
      case 4:
        GetLittleEndianFlag();
        break;
      case 10:
        GetNodesAscii();
        break;
      case 12:
        GetCellsAscii();
        break;
      case 13:
        GetFacesAscii();
        break;
      case 18:
        GetPeriodicShadowFacesAscii();
        break;
      case 37:
        break;
      case 38:
        break;
      case 39:
        break;
      case 40:
        break;
      case 41:
        break;
      case 45:
        break;
      case 58:
        GetCellTreeAscii();
        break;
      case 59:
        GetFaceTreeAscii();
        break;
      case 61:
        GetInterfaceFaceParentsAscii();
        break;
      case 62:
        GetNonconformalGridInterfaceFaceInformationAscii();
        break;
      case 63:
        break;
      case 64:
        break;
      case 2010:
        GetNodesSinglePrecision();
        break;
      case 3010:
        GetNodesDoublePrecision();
        break;
      case 2012:
        GetCellsBinary();
        break;
      case 3012:
        GetCellsBinary();  // Should be the same as single precision.. only grabbing ints.
        break;
      case 2013:
        GetFacesBinary();
        break;
      case 3013:
        GetFacesBinary();
        break;
      case 2018:
        GetPeriodicShadowFacesBinary();
        break;
      case 3018:
        GetPeriodicShadowFacesBinary();
        break;
      case 2040:
        break;
      case 3040:
        break;
      case 2041:
        break;
      case 3041:
        break;
      case 2058:
        GetCellTreeBinary();
        break;
      case 3058:
        GetCellTreeBinary();
        break;
      case 2059:
        GetFaceTreeBinary();
        break;
      case 3059:
        GetFaceTreeBinary();
        break;
      case 2061:
        GetInterfaceFaceParentsBinary();
        break;
      case 3061:
        GetInterfaceFaceParentsBinary();
        break;
      case 2062:
        GetNonconformalGridInterfaceFaceInformationBinary();
        break;
      case 3062:
        GetNonconformalGridInterfaceFaceInformationBinary();
        break;
      case 2063:
        break;
      case 3063:
        break;
      default:
        //cout << "Undefined Section = " << index << endl;
        break;
      }
    }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetDimension()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+4, 1 );
  return atoi(info.c_str());
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetLittleEndianFlag()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int flag;
  sscanf(info.c_str(), "%d", &flag);

  if (flag == 60)
    {
    LittleEndianFlag = 1;
    }
  else
    {
    LittleEndianFlag = 0;
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNodesAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int zoneId, firstIndex, lastIndex, type, nd;
  sscanf(info.c_str(), "%x %x %x %d %d", &zoneId, &firstIndex, &lastIndex, &type, &nd);

  if (CaseBuffer.at(5) == '0')
    {
    Points->Allocate(lastIndex);
    }
  else
    {
    int dstart = CaseBuffer.find('(', 5);
    int dend = CaseBuffer.find(')', dstart+1);
    vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
    vtkstd::istringstream pdatastream;
    pdatastream.str(pdata);

    double x, y, z;
    if (GridDimension == 3)
      {
      for (int i = firstIndex; i <=lastIndex; i++)
        {
        pdatastream >> x;
        pdatastream >> y;
        pdatastream >> z;
        Points->InsertPoint(i-1, x, y, z);
        }
      }
    else
      {
      for (int i = firstIndex; i <=lastIndex; i++)
        {
        pdatastream >> x;
        pdatastream >> y;
        Points->InsertPoint(i-1, x, y, 0.0);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNodesSinglePrecision()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int zoneId, firstIndex, lastIndex, type, nd;
  sscanf(info.c_str(), "%x %x %x %d", &zoneId, &firstIndex, &lastIndex, &type);

  int dstart = CaseBuffer.find('(', 7);
  int ptr = dstart + 1;

  double x, y, z;
  if (GridDimension == 3)
    {
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      x = GetCaseBufferFloat(ptr);
      ptr = ptr + 4;

      y = GetCaseBufferFloat(ptr);
      ptr = ptr + 4;

      z = GetCaseBufferFloat(ptr);
      ptr = ptr + 4;
      Points->InsertPoint(i-1, x, y, z);
      }
    }
  else
    {
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      x = GetCaseBufferFloat(ptr);
      ptr = ptr + 4;

      y = GetCaseBufferFloat(ptr);
      ptr = ptr + 4;

      z = 0.0;

      Points->InsertPoint(i-1, x, y, 0.0);
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNodesDoublePrecision()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int zoneId, firstIndex, lastIndex, type, nd;
  sscanf(info.c_str(), "%x %x %x %d", &zoneId, &firstIndex, &lastIndex, &type);

  int dstart = CaseBuffer.find('(', 7);
  int ptr = dstart+1;

  double x, y, z;
  if (GridDimension == 3)
    {
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      x = GetCaseBufferDouble(ptr);
      ptr = ptr + 8;

      y = GetCaseBufferDouble(ptr);
      ptr = ptr + 8;

      z = GetCaseBufferDouble(ptr);
      ptr = ptr + 8;
      Points->InsertPoint(i-1, x, y, z);
      }
    }
  else
    {
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      x = GetCaseBufferDouble(ptr);
      ptr = ptr + 8;

      y = GetCaseBufferDouble(ptr);
      ptr = ptr + 8;

      z = 0.0;
      Points->InsertPoint(i-1, x, y, 0.0);
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetCellsAscii()
{
  if (CaseBuffer.at(5) == '0')
    { // Cell Info
    int start = CaseBuffer.find('(', 1);
    int end = CaseBuffer.find(')',1);
    vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
    int zoneId, firstIndex, lastIndex, type;
    sscanf(info.c_str(), "%x %x %x %d", &zoneId, &firstIndex, &lastIndex, &type);
    Cells.resize(lastIndex);
    }
  else
    { // Cell Definitions
    int start = CaseBuffer.find('(', 1);
    int end = CaseBuffer.find(')',1);
    vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
    int zoneId, firstIndex, lastIndex, type, elementType;
    sscanf(info.c_str(), "%x %x %x %d %d", &zoneId, &firstIndex, &lastIndex, &type, &elementType);

    if (elementType == 0)
      {
      int dstart = CaseBuffer.find('(', 5);
      int dend = CaseBuffer.find(')', dstart+1);
      vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
      vtkstd::istringstream pdatastream;
      pdatastream.str(pdata);
      for (int i = firstIndex; i <=lastIndex; i++)
        {
        pdatastream >> Cells[i].type;
        Cells[i-1].zone = zoneId;
        Cells[i-1].parent = 0;
        Cells[i-1].child  = 0;
        }
      }
    else 
      {
      for (int i = firstIndex; i <=lastIndex; i++)
        {
        Cells[i-1].type = elementType;
        Cells[i-1].zone = zoneId;
        Cells[i-1].parent = 0;
        Cells[i-1].child  = 0;
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetCellsBinary()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int zoneId, firstIndex, lastIndex, type, elementType;
  sscanf(info.c_str(), "%x %x %x %x %x", &zoneId, &firstIndex, &lastIndex, &type, &elementType);

  if (elementType == 0)
    {
    int dstart = CaseBuffer.find('(', 7);
    int ptr = dstart + 1;
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      Cells[i-1].type = GetCaseBufferInt(ptr);
      ptr = ptr +4;
      Cells[i-1].zone = zoneId;
      Cells[i-1].parent = 0;
      Cells[i-1].child  = 0;
      }
    }
  else 
    {
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      Cells[i-1].type = elementType;
      Cells[i-1].zone = zoneId;
      Cells[i-1].parent = 0;
      Cells[i-1].child  = 0;
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetFacesAscii()
{

  if (CaseBuffer.at(5) == '0')
    { // Face Info
    int start = CaseBuffer.find('(', 1);
    int end = CaseBuffer.find(')',1);
    vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
    int zoneId, firstIndex, lastIndex, bcType, faceType;
    sscanf(info.c_str(), "%x %x %x %x", &zoneId, &firstIndex, &lastIndex, &bcType);

    Faces.resize(lastIndex);
    }
  else
    { // Face Definitions
    int start = CaseBuffer.find('(', 1);
    int end = CaseBuffer.find(')',1);
    vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
    int zoneId, firstIndex, lastIndex, bcType, faceType;
    sscanf(info.c_str(), "%x %x %x %x %x", &zoneId, &firstIndex, &lastIndex, &bcType, &faceType);

    int dstart = CaseBuffer.find('(', 7);
    int dend = CaseBuffer.find(')', dstart+1);
    vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
    vtkstd::istringstream pdatastream;
    pdatastream.str(pdata);

    int numberOfNodesInFace = 0;
    for (int i = firstIndex; i <=lastIndex; i++)
      {
      if (faceType == 0 || faceType == 5)
        {
        pdatastream >> numberOfNodesInFace;
        }
      else
        {
        numberOfNodesInFace = faceType;
        }
      Faces[i-1].nodes.resize(numberOfNodesInFace);
      for (int j = 0; j<numberOfNodesInFace; j++)
        {
        pdatastream >> hex >> Faces[i-1].nodes[j];
        Faces[i-1].nodes[j]--;
        }
      pdatastream >> hex >> Faces[i-1].c0;
      pdatastream >> hex >> Faces[i-1].c1;
      Faces[i-1].c0--;
      Faces[i-1].c1--;
      Faces[i-1].type = numberOfNodesInFace;
      Faces[i-1].zone = zoneId;
      Faces[i-1].periodicShadow = 0;
      Faces[i-1].parent = 0;
      Faces[i-1].child = 0;
      Faces[i-1].interfaceFaceParent = 0;
      Faces[i-1].ncgParent = 0;
      Faces[i-1].ncgChild = 0;
      Faces[i-1].interfaceFaceChild = 0;
      if (Faces[i-1].c0 >= 0)
        {
        Cells[Faces[i-1].c0].faces.push_back(i-1);
        }
      if (Faces[i-1].c1 >= 0)
        {
        Cells[Faces[i-1].c1].faces.push_back(i-1);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetFacesBinary()
{
//cout<<"Get Faces Binary"<<endl;
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int zoneId, firstIndex, lastIndex, bcType, faceType;
  sscanf(info.c_str(), "%x %x %x %x %x", &zoneId, &firstIndex, &lastIndex, &bcType, &faceType);
//cout<<zoneId<<" : "<<firstIndex<<" : "<<lastIndex<<" : "<<bcType<<" : "<< faceType  <<endl;
  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int numberOfNodesInFace = 0;
  int ptr = dstart + 1;
//cout<<lastIndex<<endl;
  for (int i = firstIndex; i <=lastIndex; i++)
    {
//cout<<i<<endl;
    if ((faceType == 0) || (faceType == 5))
      {
      numberOfNodesInFace = GetCaseBufferInt(ptr);
      ptr = ptr + 4;
      }
    else
      {
      numberOfNodesInFace = faceType;
      }

    Faces[i-1].nodes.resize(numberOfNodesInFace);

    for (int k = 0; k<numberOfNodesInFace; k++)
      {
      Faces[i-1].nodes[k] = GetCaseBufferInt(ptr);
      Faces[i-1].nodes[k]--;
      ptr = ptr + 4;
      }

    Faces[i-1].c0 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    Faces[i-1].c1 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    Faces[i-1].c0--;
    Faces[i-1].c1--;
    Faces[i-1].type = numberOfNodesInFace;
    Faces[i-1].zone = zoneId;
    Faces[i-1].periodicShadow = 0;
    Faces[i-1].parent = 0;
    Faces[i-1].child = 0;
    Faces[i-1].interfaceFaceParent = 0;
    Faces[i-1].ncgParent = 0;
    Faces[i-1].ncgChild = 0;
    Faces[i-1].interfaceFaceChild = 0;
//cout<<Faces[i-1].c0<<" : "<<Faces[i-1].c1<<endl;
//cout<<Cells.size()<<endl;
    if (Faces[i-1].c0 >= 0)
      {
      Cells[Faces[i-1].c0].faces.push_back(i-1);
      }
    if (Faces[i-1].c1 >= 0)
      {
      Cells[Faces[i-1].c1].faces.push_back(i-1);
      }
    }
//cout<<"end"<<endl;
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetPeriodicShadowFacesAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int firstIndex, lastIndex, periodicZone, shadowZone;
  sscanf(info.c_str(), "%x %x %x %x", &firstIndex, &lastIndex, &periodicZone, &shadowZone);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
  vtkstd::istringstream pdatastream;
  pdatastream.str(pdata);

  int faceIndex1, faceIndex2;
  for (int i = firstIndex; i <=lastIndex; i++)
    {
    pdatastream >> hex >> faceIndex1;
    pdatastream >> hex >> faceIndex2;
    Faces[faceIndex1].periodicShadow = 1;
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetPeriodicShadowFacesBinary()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int firstIndex, lastIndex, periodicZone, shadowZone;
  sscanf(info.c_str(), "%x %x %x %x", &firstIndex, &lastIndex, &periodicZone, &shadowZone);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int ptr = dstart + 1;

  int faceIndex1, faceIndex2;
  for (int i = firstIndex; i <=lastIndex; i++)
    {
    faceIndex1 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    faceIndex2 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetCellTreeAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int cellId0, cellId1, parentZoneId, childZoneId;
  sscanf(info.c_str(), "%x %x %x %x", &cellId0, &cellId1, &parentZoneId, &childZoneId);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
  vtkstd::istringstream pdatastream;
  pdatastream.str(pdata);

  int numberOfKids, kid;
  for (int i = cellId0; i <=cellId1; i++)
    {
    Cells[i-1].parent = 1;
    pdatastream >> hex >> numberOfKids;
    for (int j = 0; j < numberOfKids; j++)
      {
      pdatastream >> hex >> kid;
      Cells[kid-1].child = 1;
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetCellTreeBinary()
{

  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int cellId0, cellId1, parentZoneId, childZoneId;
  sscanf(info.c_str(), "%x %x %x %x", &cellId0, &cellId1, &parentZoneId, &childZoneId);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int ptr = dstart + 1;

  int numberOfKids, kid;
  for (int i = cellId0; i <=cellId1; i++)
    {
    Cells[i-1].parent = 1;
    numberOfKids = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    for (int j = 0; j < numberOfKids; j++)
      {
      kid = GetCaseBufferInt(ptr);
      ptr = ptr + 4;
      Cells[kid-1].child = 1;
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetFaceTreeAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int faceId0, faceId1, parentZoneId, childZoneId;
  sscanf(info.c_str(), "%x %x %x %x", &faceId0, &faceId1, &parentZoneId, &childZoneId);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
  vtkstd::istringstream pdatastream;
  pdatastream.str(pdata);

  int numberOfKids, kid;
  for (int i = faceId0; i <=faceId1; i++)
    {
    Faces[i-1].parent = 1;
    pdatastream >> hex >> numberOfKids;
    for (int j = 0; j < numberOfKids; j++)
      {
      pdatastream >> hex >> kid;
      Faces[kid-1].child = 1;
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetFaceTreeBinary()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int faceId0, faceId1, parentZoneId, childZoneId;
  sscanf(info.c_str(), "%x %x %x %x", &faceId0, &faceId1, &parentZoneId, &childZoneId);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int ptr = dstart + 1;

  int numberOfKids, kid;
  for (int i = faceId0; i <=faceId1; i++)
    {
    Faces[i-1].parent = 1;
    numberOfKids = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    for (int j = 0; j < numberOfKids; j++)
      {
      kid = GetCaseBufferInt(ptr);
      ptr = ptr + 4;
      Faces[kid-1].child = 1;
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetInterfaceFaceParentsAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int faceId0, faceId1;
  sscanf(info.c_str(), "%x %x", &faceId0, &faceId1);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
  vtkstd::istringstream pdatastream;
  pdatastream.str(pdata);

  int parentId0, parentId1;
  for (int i = faceId0; i <=faceId1; i++)
    {
    pdatastream >> hex >> parentId0;
    pdatastream >> hex >> parentId1;
    Faces[parentId0-1].interfaceFaceParent = 1;
    Faces[parentId1-1].interfaceFaceParent = 1;
    Faces[i-1].interfaceFaceChild = 1;
    }

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetInterfaceFaceParentsBinary()
{

  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int faceId0, faceId1;
  sscanf(info.c_str(), "%x %x", &faceId0, &faceId1);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int ptr = dstart + 1;

  int parentId0, parentId1;
  for (int i = faceId0; i <=faceId1; i++)
    {
    parentId0 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    parentId1 = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    Faces[parentId0-1].interfaceFaceParent = 1;
    Faces[parentId1-1].interfaceFaceParent = 1;
    Faces[i-1].interfaceFaceChild = 1;
    }

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNonconformalGridInterfaceFaceInformationAscii()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int KidId, ParentId, NumberOfFaces;
  sscanf(info.c_str(), "%d %d %d", &KidId, &ParentId, &NumberOfFaces);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  vtkstd::string pdata = CaseBuffer.substr(dstart+1, dend-start-1);
  vtkstd::istringstream pdatastream;
  pdatastream.str(pdata);

  int child, parent;
  for (int i = 0; i < NumberOfFaces; i++)
    {
    pdatastream >> hex >> child;
    pdatastream >> hex >> parent;
    Faces[child-1].ncgChild = 1;
    Faces[parent-1].ncgParent = 1;
    }

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::GetNonconformalGridInterfaceFaceInformationBinary()
{
  int start = CaseBuffer.find('(', 1);
  int end = CaseBuffer.find(')',1);
  vtkstd::string info = CaseBuffer.substr(start+1,end-start-1 );
  int KidId, ParentId, NumberOfFaces;
  sscanf(info.c_str(), "%d %d %d", &KidId, &ParentId, &NumberOfFaces);

  int dstart = CaseBuffer.find('(', 7);
  int dend = CaseBuffer.find(')', dstart+1);
  int ptr = dstart + 1;

  int child, parent;
  for (int i = 0; i < NumberOfFaces; i++)
    {
    child = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    parent = GetCaseBufferInt(ptr);
    ptr = ptr + 4;
    Faces[child-1].ncgChild = 1;
    Faces[parent-1].ncgParent = 1;
    }

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::CleanCells()
{

vtkstd::vector<int> t;

  for (int i = 0; i < Cells.size(); i++)
    {

    if ( ((Cells[i].type == 1)  && (Cells[i].faces.size() != 3)) ||
         ((Cells[i].type == 2)  && (Cells[i].faces.size() != 4)) ||
         ((Cells[i].type == 3)  && (Cells[i].faces.size() != 4)) ||
         ((Cells[i].type == 4)  && (Cells[i].faces.size() != 6)) ||
         ((Cells[i].type == 5)  && (Cells[i].faces.size() != 5)) ||
         ((Cells[i].type == 6)  && (Cells[i].faces.size() != 5)) )
      {

      // Copy faces
      t.clear();
      for (int j = 0; j < Cells[i].faces.size(); j++)
        {
        t.push_back(Cells[i].faces[j]);
        }

      // Clear Faces
      Cells[i].faces.clear();

      // Copy the faces that are not flagged back into the cell
      for (int j = 0; j < t.size(); j++)
        {
        if ( (Faces[t[j]].child == 0 ) &&
             (Faces[t[j]].ncgChild == 0 ) &&
             (Faces[t[j]].interfaceFaceChild == 0 ))
          {
          Cells[i].faces.push_back(t[j]);
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateCellNodes()
{
  for (int i = 0; i < Cells.size(); i++)
    {
cout<<"Pop: "<<Cells[i].type<<endl;
    switch (Cells[i].type)
      {
      case 1:  // Triangle
            PopulateTriangleCell(i);
            break;

      case 2:  // Tetrahedron
            PopulateTetraCell(i);
            break;

      case 3:  // Quadrilateral
            PopulateQuadCell(i);
            break;

      case 4:  // Hexahedral
            PopulateHexahedronCell(i);
            break;

      case 5:  // Pyramid
            PopulatePyramidCell(i);
            break;

      case 6:  // Wedge
            PopulateWedgeCell(i);
            break;

      case 7:  // Polyhedron
            PopulatePolyhedronCell(i);
            break;
      }
    }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetCaseBufferInt(int ptr)
{
  union mix_i
    {
    int i;
    char c[4];
    } mi = {1};

  for (int j = 0; j < 4; j++)
    {
    if (!LittleEndianFlag)
      {
      mi.c[3 - j] = CaseBuffer.at(ptr+j);
      }
    else
      {
      mi.c[j] = CaseBuffer.at(ptr+j);
      }
    }
  return mi.i;
}

//----------------------------------------------------------------------------
float vtkFLUENTReader::GetCaseBufferFloat(int ptr)
{
  union mix_f
    {
    float f;
    char c[4];
    } mf = {1.0};

  for (int j = 0; j < 4; j++)
    {
    if (!LittleEndianFlag)
      {
      mf.c[3 - j] = CaseBuffer.at(ptr+j);
      }
    else
      {
      mf.c[j] = CaseBuffer.at(ptr+j);
      }
    }
  return mf.f;
}

//----------------------------------------------------------------------------
double vtkFLUENTReader::GetCaseBufferDouble(int ptr)
{
  union mix_i
    {
    double d;
    char c[8];
    } md = {1.0};

  for (int j = 0; j < 8; j++)
    {
    if (!LittleEndianFlag)
      {
      md.c[7 - j] = CaseBuffer.at(ptr+j);
      }
    else
      {
      md.c[j] = CaseBuffer.at(ptr+j);
      }
    }
  return md.d;
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateTriangleCell(int i)
{
  Cells[i].nodes.resize(3);
  if (Faces[Cells[i].faces[0]].c0 == i)
    {
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[1];
    }
  else
    {
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[1];
    }

  if (Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[0] &&
      Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[1])
    {
    Cells[i].nodes[2] = Faces[Cells[i].faces[1]].nodes[0];
    }
  else
    {
    Cells[i].nodes[2] = Faces[Cells[i].faces[1]].nodes[1];
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateTetraCell(int i)
{
  Cells[i].nodes.resize(4);

  if (Faces[Cells[i].faces[0]].c0 == i)
    {
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[1];
    Cells[i].nodes[2] = Faces[Cells[i].faces[0]].nodes[2];
    }
  else
    {
    Cells[i].nodes[2] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[1];
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[2];
    }

  if (Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[0] &&
      Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[1] &&
      Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[2] )
    {
    Cells[i].nodes[3] = Faces[Cells[i].faces[1]].nodes[0];
    }
  else if (Faces[Cells[i].faces[1]].nodes[1] != Cells[i].nodes[0] &&
           Faces[Cells[i].faces[1]].nodes[1] != Cells[i].nodes[1] &&
           Faces[Cells[i].faces[1]].nodes[1] != Cells[i].nodes[2] )
    {
    Cells[i].nodes[3] = Faces[Cells[i].faces[1]].nodes[1];
    }
  else
    {
    Cells[i].nodes[3] = Faces[Cells[i].faces[1]].nodes[2];
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateQuadCell(int i)
{
  Cells[i].nodes.resize(4);

  if (Faces[Cells[i].faces[0]].c0 == i)
    {
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[1];
    }
  else
    {
    Cells[i].nodes[1] = Faces[Cells[i].faces[0]].nodes[0];
    Cells[i].nodes[0] = Faces[Cells[i].faces[0]].nodes[1];
    }

  if ((Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[0] &&
       Faces[Cells[i].faces[1]].nodes[0] != Cells[i].nodes[1]) && 
      (Faces[Cells[i].faces[1]].nodes[1] != Cells[i].nodes[0] &&
       Faces[Cells[i].faces[1]].nodes[1] != Cells[i].nodes[1]))
    {
    if (Faces[Cells[i].faces[1]].c0 == i)
      {
      Cells[i].nodes[2] = Faces[Cells[i].faces[1]].nodes[0];
      Cells[i].nodes[3] = Faces[Cells[i].faces[1]].nodes[1];
      }
    else
      {
      Cells[i].nodes[3] = Faces[Cells[i].faces[1]].nodes[0];
      Cells[i].nodes[2] = Faces[Cells[i].faces[1]].nodes[1];
      }
    }
  else if ((Faces[Cells[i].faces[2]].nodes[0] != Cells[i].nodes[0] &&
            Faces[Cells[i].faces[2]].nodes[0] != Cells[i].nodes[1]) && 
           (Faces[Cells[i].faces[2]].nodes[1] != Cells[i].nodes[0] &&
            Faces[Cells[i].faces[2]].nodes[1] != Cells[i].nodes[1]))
    {
    if (Faces[Cells[i].faces[2]].c0 == i)
      {
      Cells[i].nodes[2] = Faces[Cells[i].faces[2]].nodes[0];
      Cells[i].nodes[3] = Faces[Cells[i].faces[2]].nodes[1];
      }
    else
      {
      Cells[i].nodes[3] = Faces[Cells[i].faces[2]].nodes[0];
      Cells[i].nodes[2] = Faces[Cells[i].faces[2]].nodes[1];
      }
    }
  else
    {
    if (Faces[Cells[i].faces[3]].c0 == i)
      {
      Cells[i].nodes[2] = Faces[Cells[i].faces[3]].nodes[0];
      Cells[i].nodes[3] = Faces[Cells[i].faces[3]].nodes[1];
      }
    else
      {
      Cells[i].nodes[3] = Faces[Cells[i].faces[3]].nodes[0];
      Cells[i].nodes[2] = Faces[Cells[i].faces[3]].nodes[1];
      }
    }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateHexahedronCell(int i)
{
  Cells[i].nodes.resize(8);

  if (Faces[Cells[i].faces[0]].c0 == i)
    {
    for (int j = 0; j < 4; j++)
      {
      Cells[i].nodes[j] = Faces[Cells[i].faces[0]].nodes[j];
      }
    }
  else
    {
    for (int j = 3; j >=0; j--)
      {
      Cells[i].nodes[3-j] = Faces[Cells[i].faces[0]].nodes[j];
      }
    }

  //  Look for opposite face of hexahedron
  for (int j = 1; j < 6; j++)
    {
    int flag = 0;
    for (int k = 0; k < 4; k++)
      {
      if ( (Cells[i].nodes[0] == Faces[Cells[i].faces[j]].nodes[k]) ||
           (Cells[i].nodes[1] == Faces[Cells[i].faces[j]].nodes[k]) ||
           (Cells[i].nodes[2] == Faces[Cells[i].faces[j]].nodes[k]) ||
           (Cells[i].nodes[3] == Faces[Cells[i].faces[j]].nodes[k]) )
        {
        flag = 1;
        }
      }
    if (flag == 0)
      {
      if (Faces[Cells[i].faces[j]].c1 == i)
        {
        for (int k = 4; k < 8; k++)
          {
          Cells[i].nodes[k] = Faces[Cells[i].faces[j]].nodes[k-4];
          }
        }
      else
        {
        for (int k = 7; k >= 4; k--)
          {
          Cells[i].nodes[k] = Faces[Cells[i].faces[j]].nodes[7-k];
          }
        }
      }
    }

  //  Find the face with points 0 and 1 in them.
  int f01[4];
  for (int j = 1; j < 6; j++)
    {
    int flag0 = 0;
    int flag1 = 0;
    for (int k = 0; k < 4; k++)
      {
      if (Cells[i].nodes[0] == Faces[Cells[i].faces[j]].nodes[k])
        {
        flag0 = 1;
        }
      if (Cells[i].nodes[1] == Faces[Cells[i].faces[j]].nodes[k])
        {
        flag1 = 1;
        }
      }
    if ((flag0 == 1) && (flag1 == 1))
      {
      if (Faces[Cells[i].faces[j]].c0 == i)
        {
        for (int k=0; k<4; k++)
          {
          f01[k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      else
        {
        for (int k=3; k>=0; k--)
          {
          f01[k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      }
    }

  //  Find the face with points 0 and 3 in them.
  int f03[4];
  for (int j = 1; j < 6; j++)
    {
    int flag0 = 0;
    int flag1 = 0;
    for (int k = 0; k < 4; k++)
      {
      if (Cells[i].nodes[0] == Faces[Cells[i].faces[j]].nodes[k])
        {
        flag0 = 1;
        }
      if (Cells[i].nodes[3] == Faces[Cells[i].faces[j]].nodes[k])
        {
        flag1 = 1;
        }
      }

    if ((flag0 == 1) && (flag1 == 1))
      {
      if (Faces[Cells[i].faces[j]].c0 == i)
        {
        for (int k=0; k<4; k++)
          {
          f03[k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      else
        {
        for (int k=3; k>=0; k--)
          {
          f03[k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      }
    }

  // What point is in f01 and f03 besides 0 ... this is point 4
  int p4;
  for (int k = 0; k < 4; k++)
    {
    if ( f01[k] != Cells[i].nodes[0]) 
      {
      for (int n = 0; n < 4; n++)
        {
        if (f01[k] == f03[n])
          {
          p4 = f01[k];
          }
        }
      }
    }

  // Since we know point 4 now we check to see if points
  //  4, 5, 6, and 7 are in the correct positions.
  int t[8];
  t[4] = Cells[i].nodes[4];
  t[5] = Cells[i].nodes[5];
  t[6] = Cells[i].nodes[6];
  t[7] = Cells[i].nodes[7];
  if (p4 == Cells[i].nodes[5])
    {
    Cells[i].nodes[5] = t[6];
    Cells[i].nodes[6] = t[7];
    Cells[i].nodes[7] = t[4];
    Cells[i].nodes[4] = t[5];
    }
  else if (p4 == Cells[i].nodes[6])
    {
    Cells[i].nodes[5] = t[7];
    Cells[i].nodes[6] = t[4];
    Cells[i].nodes[7] = t[5];
    Cells[i].nodes[4] = t[6];
    }
  else if (p4 == Cells[i].nodes[7])
    {
    Cells[i].nodes[5] = t[4];
    Cells[i].nodes[6] = t[5];
    Cells[i].nodes[7] = t[6];
    Cells[i].nodes[4] = t[7];
    }
  // else point 4 was lined up so everything was correct.
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulatePyramidCell(int i)
{
  Cells[i].nodes.resize(5);
  //  The quad face will be the base of the pyramid
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if ( Faces[Cells[i].faces[j]].nodes.size() == 4)
      {
      if (Faces[Cells[i].faces[j]].c0 == i)
        {
        for (int k = 0; k < 4; k++)
          {
          Cells[i].nodes[k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      else
        {
        for (int k = 0; k < 4; k++)
          {
          Cells[i].nodes[3-k] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      }
    }

  // Just need to find point 4
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if ( Faces[Cells[i].faces[j]].nodes.size() == 3)
      {
      for (int k = 0; k < 3; k ++)
        {
        if ( (Faces[Cells[i].faces[j]].nodes[k] != Cells[i].nodes[0]) &&
             (Faces[Cells[i].faces[j]].nodes[k] != Cells[i].nodes[1]) &&
             (Faces[Cells[i].faces[j]].nodes[k] != Cells[i].nodes[2]) && 
             (Faces[Cells[i].faces[j]].nodes[k] != Cells[i].nodes[3]) )
          {
          Cells[i].nodes[4] = Faces[Cells[i].faces[j]].nodes[k];
          }
        }
      }
    }

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulateWedgeCell(int i)
{
  Cells[i].nodes.resize(6);

  //  Find the first triangle face and make it the base.
  int base = 0;
  int first = 0;
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if ((Faces[Cells[i].faces[j]].type == 3) && (first == 0))
      {
      base = Cells[i].faces[j];
      first = 1;
      }
    }

  //  Find the second triangle face and make it the top.
  int top = 0;
  int second = 0;
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if ((Faces[Cells[i].faces[j]].type == 3) && (second == 0) && (Cells[i].faces[j] != base))
      {
      top = Cells[i].faces[j];
      second = 1;
      }
    }

  // Load Base nodes into the nodes vtkstd::vector
  if (Faces[base].c0 == i)
    {
    for (int j = 0; j < 3; j++)
      {
      Cells[i].nodes[j] = Faces[base].nodes[j];
      }
    }
  else
    {
    for (int j = 2; j >=0; j--)
      {
      Cells[i].nodes[2-j] = Faces[base].nodes[j];
      }
    }
   // Load Top nodes into the nodes vtkstd::vector
  if (Faces[top].c1 == i)
    {
    for (int j = 3; j < 6; j++)
      {
      Cells[i].nodes[j] = Faces[top].nodes[j-3];
      }
    }
  else
    {
    for (int j = 3; j < 6; j++)
      {
      Cells[i].nodes[j] = Faces[top].nodes[5-j];
      }
    }

  //  Find the quad face with points 0 and 1 in them.
  int w01[4];
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if (Cells[i].faces[j] != base)
      {
      int wf0 = 0;
      int wf1 = 0;
      for (int k = 0; k < 4; k++)
        {
        if (Cells[i].nodes[0] == Faces[Cells[i].faces[j]].nodes[k])
          {
          wf0 = 1;
          }
        if (Cells[i].nodes[1] == Faces[Cells[i].faces[j]].nodes[k])
          {
          wf1 = 1;
          }
        if ((wf0 == 1) && (wf1 == 1))
          {
          for (int n=0; n<4; n++)
            {
            w01[n] = Faces[Cells[i].faces[j]].nodes[n];
            }
          }
        }
      }
    }

  //  Find the quad face with points 0 and 2 in them.
  int w02[4];
  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    if (Cells[i].faces[j] != base)
      {
      int wf0 = 0;
      int wf2 = 0;
      for (int k = 0; k < 4; k++)
        {
        if (Cells[i].nodes[0] == Faces[Cells[i].faces[j]].nodes[k])
          {
          wf0 = 1;
          }
        if (Cells[i].nodes[2] == Faces[Cells[i].faces[j]].nodes[k])
          {
          wf2 = 1;
          }
        if ((wf0 == 1) && (wf2 == 1))
          {
          for (int n=0; n<4; n++)
            {
            w02[n] = Faces[Cells[i].faces[j]].nodes[n];
            }
          }
        }
      }
    }

  // Point 3 is the point that is in both w01 and w02

  // What point is in f01 and f02 besides 0 ... this is point 3
  int p3;
  for (int k = 0; k < 4; k++)
    {
    if ( w01[k] != Cells[i].nodes[0]) 
      {
      for (int n = 0; n < 4; n++)
        {
        if (w01[k] == w02[n])
          {
          p3 = w01[k];
          }
        }
      }
    }

  // Since we know point 3 now we check to see if points
  //  3, 4, and 5 are in the correct positions.
  int t[6];
  t[3] = Cells[i].nodes[3];
  t[4] = Cells[i].nodes[4];
  t[5] = Cells[i].nodes[5];
  if (p3 == Cells[i].nodes[4])
    {
    Cells[i].nodes[3] = t[4];
    Cells[i].nodes[4] = t[5];
    Cells[i].nodes[5] = t[3];
    }
  else if (p3 == Cells[i].nodes[5])
    {
    Cells[i].nodes[3] = t[5];
    Cells[i].nodes[4] = t[3];
    Cells[i].nodes[5] = t[4];
    }
  // else point 3 was lined up so everything was correct.

}

//----------------------------------------------------------------------------
void vtkFLUENTReader::PopulatePolyhedronCell(int i)
{
  //  We can't set the size on the nodes vtkstd::vector because we
  //  are not sure how many we are going to have.
  //  All we have to do here is add the nodes from the faces into
  //  nodes vtkstd::vector within the cell.  All we have to check for is 
  //  duplicate nodes.
  //
  //cout << "number of faces in cell = " << Cells[i].faces.size() << endl;

  for (int j = 0; j < Cells[i].faces.size(); j++)
    {
    //cout << "number of nodes in face = " << Faces[Cells[i].faces[j]].nodes.size() << endl;
    for (int k = 0; k < Faces[Cells[i].faces[j]].nodes.size(); k++)
      {
      int flag;
      flag = 0;
      // Is the node already in the cell?
      for (int n = 0; n < Cells[i].nodes.size(); n++)
        {
        if (Cells[i].nodes[n] == Faces[Cells[i].faces[j]].nodes[k])
          {
          flag = 1;
          }
        }
      if (flag == 0)
       { // No match - insert node into cell.
       Cells[i].nodes.push_back(Faces[Cells[i].faces[j]].nodes[k]);
      // cout << "insert node" << endl;
       }
     }
   }
}

//----------------------------------------------------------------------------
void vtkFLUENTReader::ParseDataFile()
{
  int cnt = 0;
  while (GetDataChunk())
    {
    int index = GetDataIndex();
cout<<"parsedata: "<<index<<endl;
    switch (index)
      {
      case 0:
        //cout << "Comment Section" << endl;
        break;

      case 4:
        //cout << "Machine Configuration Section" << endl;
        break;

      case 33:
        //cout << "Grid Size Section" << endl;
        break;

      case 37:
        //cout << "Variables Section" << endl;
        break;

      case 300:
        //cout << "Data Section" << endl;
        GetData(1);
        break;

      case 301:
        //cout << "Residuals Section" << endl;
        break;

      case 302:
        //cout << "Residuals Section" << endl;
        break;

      case 2300:
        //cout << "Single Precision Data Section" << endl;
        GetData(2);
        break;

      case 2301:
        //cout << "Single Precision Residuals Section" << endl;
        break;

      case 2302:
        //cout << "Single Precision Residuals Section" << endl;
        break;

      case 3300:
        //cout << "Single Precision Data Section" << endl;
        GetData(3);
        break;

      case 3301:
        //cout << "Single Precision Residuals Section" << endl;
        break;

      case 3302:
        //cout << "Single Precision Residuals Section" << endl;
        break;

      default:
        //cout << "Data Undefined Section = " << index << endl;
        break;
      }
    }
}

//----------------------------------------------------------------------------
int vtkFLUENTReader::GetDataBufferInt(int ptr)
{
  union mix_i
    {
    int i;
    char c[4];
    } mi = {1};

  for (int j = 0; j < 4; j++)
    {
    if (!LittleEndianFlag)
      {
      mi.c[3 - j] = DataBuffer.at(ptr+j);
      }
    else
      {
      mi.c[j] = DataBuffer.at(ptr+j);
      }
    }
  return mi.i;
}

//----------------------------------------------------------------------------
float vtkFLUENTReader::GetDataBufferFloat(int ptr)
{
  union mix_f
    {
    float f;
    char c[4];
    } mf = {1.0};

  for (int j = 0; j < 4; j++)
    {
    if (!LittleEndianFlag)
      {
      mf.c[3 - j] = DataBuffer.at(ptr+j);
      }
    else
      {
      mf.c[j] = DataBuffer.at(ptr+j);
      }
    }
  return mf.f;
}

//----------------------------------------------------------------------------
double vtkFLUENTReader::GetDataBufferDouble(int ptr)
{
  union mix_i
    {
    double d;
    char c[8];
    } md = {1.0};

  for (int j = 0; j < 8; j++)
    {
    if (!LittleEndianFlag)
      {
      md.c[7 - j] = DataBuffer.at(ptr+j);
      }
    else
      {
      md.c[j] = DataBuffer.at(ptr+j);
      }
    }
  return md.d;
}

//------------------------------------------------------------------------------
void vtkFLUENTReader::GetData(int dataType)
{
  int start = DataBuffer.find('(', 1);
  int end = DataBuffer.find(')',1);
  vtkstd::string info = DataBuffer.substr(start+1,end-start-1 );
  vtkstd::istringstream infostream;
  infostream.str(info);
  int subSectionId, zoneId, size, nTimeLevels, nPhases, firstId, lastId;
  infostream >> subSectionId >> zoneId >> size >> nTimeLevels >> nPhases >> firstId >> lastId;

  // Is this a cell zone?
  int zmatch = 0;
  for (int i = 0; i < CellZones.size(); i++)
    {
    if (CellZones[i] == zoneId)
      {
      zmatch = 1;
      }
    }

  if (zmatch)
    {

    // Set up stream or pointer to data
    int dstart = DataBuffer.find('(', 7);
    int dend = DataBuffer.find(')', dstart+1);
    vtkstd::istringstream pdatastream;
    vtkstd::string pdata = DataBuffer.substr(dstart+1, dend-dstart-2);
    pdatastream.str(pdata);
    int ptr = dstart + 1;

    // Is this a new variable?
    int match = 0;
    for (int i = 0; i < SubSectionIds.size(); i++)
      {
      if (subSectionId == SubSectionIds[i])
        {
        match = 1;
        }
      }

    if ((match == 0) && (size < 4))
      { // new variable
      SubSectionIds.push_back(subSectionId);
      SubSectionSize.push_back(size);
      SubSectionZones.resize(SubSectionZones.size()+1);
      SubSectionZones[SubSectionZones.size()-1].push_back(zoneId);
      }

    if (size == 1)
      {
      NumberOfScalars++;
      ScalarDataChunks.resize(ScalarDataChunks.size() + 1);
      ScalarDataChunks[ScalarDataChunks.size()-1].subsectionId = subSectionId;
      ScalarDataChunks[ScalarDataChunks.size()-1].zoneId = zoneId;
      for (int i=firstId; i<=lastId; i++)
        {
        double temp;
        if (dataType == 1)
          {
          pdatastream >> temp;
          }
        else if (dataType == 2)
          {
          temp = GetDataBufferFloat(ptr);
          ptr = ptr + 4;
          }
        else
          {
          temp = GetDataBufferDouble(ptr);
          ptr = ptr + 8;
          }
        ScalarDataChunks[ScalarDataChunks.size()-1].scalarData.push_back(temp);
        }
      }
    else if (size == 3)
      {
      NumberOfVectors++;
      VectorDataChunks.resize(VectorDataChunks.size() + 1);
      VectorDataChunks[VectorDataChunks.size() - 1].subsectionId = subSectionId;
      VectorDataChunks[VectorDataChunks.size() - 1].zoneId = zoneId;
      for (int i=firstId; i<=lastId; i++)
        {
        double tempx, tempy, tempz;

        if (dataType == 1)
          {
          pdatastream >> tempx;
          pdatastream >> tempy;
          pdatastream >> tempz;
          }
        else if (dataType == 2)
          {
          tempx = GetDataBufferFloat(ptr);
          ptr = ptr + 4;
          tempy = GetDataBufferFloat(ptr);
          ptr = ptr + 4;
          tempz = GetDataBufferFloat(ptr);
          ptr = ptr + 4;
          }
        else
          {
          tempx = GetDataBufferDouble(ptr);
          ptr = ptr + 8;
          tempy = GetDataBufferDouble(ptr);
          ptr = ptr + 8;
          tempz = GetDataBufferDouble(ptr);
          ptr = ptr + 8;
          }
        VectorDataChunks[VectorDataChunks.size()-1].iComponentData.push_back(tempx);
        VectorDataChunks[VectorDataChunks.size()-1].jComponentData.push_back(tempy);
        VectorDataChunks[VectorDataChunks.size()-1].kComponentData.push_back(tempz);
        }
      }
    else
      {
      //cout << "Weird Variable Size = " << size << endl;
      }
    }
}
