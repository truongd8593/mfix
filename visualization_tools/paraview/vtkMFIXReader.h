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
// .NAME vtkMFIXReader - reads a dataset in MFIX file format
// .SECTION Description
// vtkMFIXReader creates an unstructured grid dataset. It reads a restart
// file and a set of sp files.  The restart file contains the mesh 
// information.  MFIX meshes are either cylindrical or rectilinear, but 
// this reader will convert them to an unstructured grid.  The sp files 
// contain transient data for the cells.  Each sp file has one or more 
// variables stored inside it.  

// .SECTION Thanks
// Thanks to Phil Nicoletti and Brian Dotson at the National Energy 
// Technology Laboratory who developed this class.
// Please address all comments to Phil Nicoletti (philip.nicoletti@pp.netl.doe.gov)
// or Brian Dotson (brian.dotson@netl.doe.gov)

// .SECTION See Also
// vtkGAMBITReader

#ifndef __vtkMFIXReader_h
#define __vtkMFIXReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkDataArraySelection.h"
#include "vtkStringArray.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkWedge.h"
#include "vtkHexahedron.h"

class vtkDoubleArray;
class VTK_IO_EXPORT vtkMFIXReader : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkMFIXReader *New();
  vtkTypeRevisionMacro(vtkMFIXReader,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Methods for Paraview
  int GetNumberOfCellArrays(void);
  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();
  void GetCellDataRange(int cellComp, int index, float *min, float *max);
  
  // Description:
  // Specify the file name of the MFIX Restart data file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the total number of cells. The number of cells is only valid after a
  // successful read of the data file is performed.
  vtkGetMacro(NumberOfCells,int);

  // Description:
  // Get the total number of nodes. The number of nodes is only valid after a
  // successful read of the data file is performed.
  vtkGetMacro(NumberOfPoints,int);

  // Description:
  // Get the number of data components at the nodes and cells.
  vtkGetMacro(NumberOfCellFields,int);
  
  // Description:
  // Which TimeStep to read.    
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
 
  vtkGetMacro(NumberOfTimeSteps, int);
  // Description:
  // Which TimeStepRange to read
  vtkGetVector2Macro(TimeStepRange, int);
  vtkSetVector2Macro(TimeStepRange, int);

protected:
  vtkMFIXReader();
  ~vtkMFIXReader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  //
  // ParaView Variables
  //
  
  char *FileName;
  int RequestInformationFlag;
  int MakeMeshFlag;
  int NumberOfPoints;
  int NumberOfCells;
  int NumberOfCellFields;
  int *veclen;
  float *min;
  float *max;
  vtkDataArraySelection *CellDataArraySelection;
  int TimeStep;
  int CurrentTimeStep;
  int NumberOfTimeSteps;
  int *TimeSteps; 
  int TimeStepRange[2];
  int TimeStepWasReadOnce;
  
  //
  //  MFIX Variables
  //

  vtkFloatArray **cell_data_array; // Arrays for variables that will attach to mesh
  vtkPoints *points;		// Points array for building grid
  vtkUnstructuredGrid *mesh;	// Unstructured Grid
  vtkHexahedron *aHexahedron;	// Hexahedron type cell
  vtkWedge *aWedge;  		// Wedge type cell
  vtkIntArray *FLAG;		// Cell Flag array
  vtkDoubleArray *DX;		// Cell widths in x axis
  vtkDoubleArray *DY;		// Cell widths in y axis
  vtkDoubleArray *DZ;		// Cell widths in z axis
  vtkIntArray *NMAX;		// Array to hold number of species per phase
  vtkDoubleArray *C;		// Array used to parse restart file
  vtkIntArray *tmpI;		// Array used to parse restart file
  vtkDoubleArray *tmpD;		// Array used to parse restart file
  vtkIntArray *SpxFileExists;   // Array for keeping track of what spx files exist.
  
  char ext [15];
  char buffer[513];
  char version[120];
  float version_number;
  double P_ref;
  double P_scale;
  int    DIM_IC;
  int    DIM_BC;
  int    DIM_C;
  int    DIM_IS;
  double c_e;
  double c_f;
  double phi;
  double phi_w;
  double C_e;
  double C_f;
  double Phi;
  double Phi_w;
  double dt;
  double xmin;
  char run_name[256];
  vtkStringArray *variable_names;
  vtkIntArray *variable_components;
  int imin1;
  int jmin1;
  int kmin1;
  int imax;
  int jmax;
  int kmax;
  int imax1;
  int jmax1;
  int kmax1;
  int imax2;
  int jmax2;
  int kmax2;
  int ijmax2;
  int ijkmax2;
  int MMAX;
  int nspx_use;
  double xlength;
  double ylength;
  double zlength;
  int  NScalar;
  int nRR;
  int nTimes;
  bool bKepsilon;
  char coordinates[17];
  char units[17];

  //
  //  SPX Variables
  //
	
  int max_timestep;                 // maximum timesteps amongst the variables
  int spx_records_per_timestep;     // number of records in a single timestep for a variable
  vtkIntArray *spx_to_nvar_table;   // number of variables in each spx file
  vtkIntArray *var_to_skip_table;   // skip value for each variable, this is needed in spx files 
  				    // with more than one variable.
  vtkIntArray *var_timesteps;       // number of timesteps for each variable
  vtkIntArray *var_timestep_table;  //  Since the number of timesteps vary between variables
  				    //  this is a table that looks up the appropriate 
				    //  timestep for the particular variable.
  vtkIntArray *variableIndexToSPX;  //  This gives the spx file number for the 
  				    //  particular variable.
  int *spx_timestep_index_table;    //  This a table look up for the index into a file for 
  				    //  a certain variable.


private:
  void MakeMesh(vtkUnstructuredGrid *output);
  void SWAP_DOUBLE(double &value);
  void SWAP_FLOAT(float &value);
  void SWAP_INT(int &value);
  vtkStdString ConvertIntToString(int in);
  int ConvertCharToInt(char in);
  int ConvertStringToInt(const vtkStdString & in);
  void GetInt(istream& in, int &val);
  void GetDouble(istream& in, double& val);
  void GetFloat(istream& in, float& val);
  void SkipBytes(istream& in, int n);
  void RestartVersionNumber(char* buffer);
  void IN_BIN_512(istream& in, vtkDoubleArray *v, int n);
  void IN_BIN_512R(istream& in, vtkFloatArray *v, int n);
  void IN_BIN_512I(istream& in, vtkIntArray *v, int n);
  void ReadRes0();
  void GetVariableAtTimestep(int vari , int tstep, vtkFloatArray *v);
  void CreateVariableNames();
  void GetTimeSteps();
  void MakeTimeStepTable(int nvars);
  void SetProjectName (char *infile);
  void MakeSPXTimeStepIndexTable(int nvars);
  void CalculateMaxTimeStep();
  void GetNumberOfVariablesInSPXFiles();
  void FillVectorVariable( int xindex, int yindex, int zindex, vtkFloatArray *v);
   
};

#endif
