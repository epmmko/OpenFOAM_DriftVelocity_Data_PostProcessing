{
    #include "continuityErrs.H"

    wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll (p_rgh.boundaryField(), i)
    {
        if (p_rgh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

    adjustPhi(phi, U, pcorr);

    for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (nonOrth == pimple.nNonOrthCorr())
        {
            phi -= pcorrEqn.flux();
        }
    }

    #include "continuityErrs.H"
}
