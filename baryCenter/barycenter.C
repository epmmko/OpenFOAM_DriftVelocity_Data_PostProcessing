/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    barycenter

Description
    To calculate the mass center of a bubble (VOF=0) in the liquid 
    Derived from interFoam
    Changed in createFields.H: gamma:NO_WRITE

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "interpolationTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createTimeControls.H"
#   include "createMesh.H"
//#   include "readEnvironmentalProperties.H"

pimpleControl pimple(mesh);

//#   include "readPISOControls.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "readTimeControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    scalar x = 0.0;
    scalar sumM = 0.0;
    scalar sumMX = 0.0;
    //scalar oldy = 0.0;
    const scalar& write_interval(
	  readScalar(runTime.controlDict().lookup("writeInterval"))
	  );

    runTime.setDeltaT
    (
        min
        (
            write_interval,
            maxDeltaT
        )
    );


    while (runTime.run())
    {
      runTime++;
      Info << "deltaT = "<<  write_interval << endl;

      Info<< "Time = " << runTime.timeName() << nl << endl;
       
      // update the alpha1 field
      volScalarField alpha
	(
	 IOobject
	 (
	  "alpha",
	  runTime.timeName(),
	  mesh,
	  IOobject::READ_IF_PRESENT,
	  IOobject::NO_WRITE
	  ),
	 mesh
	 );

      // calc barycenter of the air bubble (phase 2, alpha1(air phase)=0!!)
      x = 0.0;
      sumM = 0.0;
      sumMX = 0.0;

      forAll(alpha, cellI){
	  if (alpha[cellI] < 2){ 
		// Integration of {mass}
	  	sumM += (1-alpha[cellI])*mesh.V()[cellI]*rho2.value();

	  	// Integration of (mass*height)
	  	sumMX += (1-alpha[cellI])*mesh.V()[cellI]*rho2.value()*mesh.C()[cellI][0];
	  }
	}

	x = sumMX/max(sumM, SMALL); // barycenter = ...

	//#include "writeBaryCenter.H"
 
	Info << "barycenter at time " <<  runTime.timeName() << " s"
	     << " is " << x << " m"<< endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
