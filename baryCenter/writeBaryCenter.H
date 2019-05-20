    Info<< "writing barycenter\n" << endl;


    volScalarField x1 = (x*rho)/rho;

	x1.write();

	volScalarField bary
	(
	 IOobject
	 (
	  "bary",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
	  ),
	 x1
	 );

	graph gr
	(
		bary_center,
		time,
		center,
		x1
	);

	gr.writeTable();
		
