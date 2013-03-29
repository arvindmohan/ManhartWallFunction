/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ManhartWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// A patch is a list of labels that address the faces in the global face list

ManhartWallFunctionFvPatchScalarField::
ManhartWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    kappa_(0.41),
    E_(9.8)
{}


ManhartWallFunctionFvPatchScalarField::
ManhartWallFunctionFvPatchScalarField
(
    const ManhartWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


ManhartWallFunctionFvPatchScalarField::
ManhartWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


ManhartWallFunctionFvPatchScalarField::
ManhartWallFunctionFvPatchScalarField
(
    const ManhartWallFunctionFvPatchScalarField& nwfpsf
)
:
    fixedValueFvPatchScalarField(nwfpsf),
    kappa_(nwfpsf.kappa_),
    E_(nwfpsf.E_)
{}


ManhartWallFunctionFvPatchScalarField::
ManhartWallFunctionFvPatchScalarField
(
    const ManhartWallFunctionFvPatchScalarField& nwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(nwfpsf, iF),
    kappa_(nwfpsf.kappa_),
    E_(nwfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ManhartWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    // Refer fvPatch.C for boundaryField,DeltaCoeffs etc.
    // Refer lesModel.H and .C
    // DeltaCoeffs: Return the face - cell distance coefficient except for coupled patches for which the cell-centre to coupled-cell-centre 
    // distance coefficient is returned
    const LESModel& lesModel = db().lookupObject<LESModel>("LESProperties"); // Read LES model information from LES properties file
    const label patchi = patch().index();  //- Return the index of this patch in the fvBoundaryMesh, and name the function as a label "patchi"
    const fvPatchVectorField& U = lesModel.U().boundaryField()[patchi]; //Define "U" as a pointer to a velocity field array as a type 
                                                          // given by "fvPatchVectorField" , and the ID of the patch is given by patchi
    const scalarField nuw = lesModel.nu()().boundaryField()[patchi]; // Same array as above, but for viscosity
    const volScalarField& pr = this->db().lookupObject<volScalarField>("p"); // Read "Pressure" from the database structure -WHATS ITS ARRAY SIZE? =  COMPLETE FIELD
    const scalarField& ry = patch().deltaCoeffs(); // Create an array having the reciprocal of the "y" distance values of boundary cells 
                                                   // from wall with LES Delta Coefficient values for all cells in the patch

    const scalarField magUp(mag(U.patchInternalField() - U)); // U(Cells Adjacent to boundary Field) - U(boundary field) --  WHY??
    scalarField& nuSgsw = *this; // Consider value of nuSgs in "this" (the current) instance. Using *this operator. nuSGS is calc. by LES model
                                 // and its value is modified in this code to suit velocity profile

    const scalarField magFaceGradU(mag(U.snGrad()));  // Calc. Surface Normal Gradient of the boundary field

    volVectorField gradp  = fvc::grad(pr);
    vectorField Gbpr = gradp.boundaryField()[patchi];	
   // const  fvPatchScalarField& bpr = pr.boundaryField()[patchi]; 
    
   //  volVectorField Gbpr = fvc::grad(bpr);
    const scalarField GradP = Gbpr.component(0); 

    forAll(nuSgsw, facei) // facei for how many faces ? = boundary field?
    {
        scalar magUpara = magUp[facei];  // Wall Parallel "U" velocity at the cell adjacent to the wall
        
        scalar tau = (nuSgsw[facei] + nuw[facei])*magFaceGradU[facei];  // Wall shear stress (tau) definition.        

        scalar utau = sqrt(mag(tau)); // Friction Velocity (Utau) definition. 

        scalar uP =  pow(mag((nuSgsw[facei] + nuw[facei])*(GradP[facei])),0.3333333333333333);  // Defining streamwise pressure based velocity (nuw or nusgs?)

        scalar utauP = sqrt(sqr(utau) + sqr(uP));

        scalar alpha = sqr(utau)/sqr(utauP);

        scalar Ustar = magUpara/utauP;

        scalar  ystar =  utauP/(ry[facei]*nuw[facei]);   

        scalar term2 =  1.0/(ry[facei]*nuw[facei]);

        if (utauP > VSMALL)   
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {  //Newton-Raphson Solution
               
               scalar f = sign(GradP[facei])*(pow((1-alpha),1.5)/2.0)*sqr(ystar) + sign(tau)*alpha*ystar - Ustar;

               scalar term = sqr(utauP) - sqr(utau);

               scalar df = sign(GradP[facei])*((3*sqr(utauP)*pow(term,0.5) - pow(term,1.5))/(2*sqr(utauP)))*sqr(term2)     
                    + (magUpara/sqr(utauP)) - sign(tau)*alpha*term2;


               scalar utauPNEW = utauP - f/df;
 
               scalar err = mag((utauP - utauPNEW)/utauP);

               utauP = utauPNEW;
             
            } while (utauP > VSMALL && err > 0.01 && ++iter < 20);        // Setting tolerance criteria for Utau and max limit of NR iterations

            // Calculating final parameters
             
             utau = sqrt(sqr(utauP) - sqr(uP));

             tau = sqr(utau); 

             scalar nuCorr = (tau/magFaceGradU[facei])*(1 + (nuSgsw[facei]/nuw[facei])) - nuw[facei] - nuSgsw[facei]; 

             nuSgsw[facei] = nuSgsw[facei] + nuCorr; // Updating new value of nuSgs

        }
        else
        {
            nuSgsw[facei] = 0;       //Seting SGS viscosity to 0 if Wall Shear stress is lesser than the VSMALL value set in OpenFOAM
        }
    }

    fixedValueFvPatchScalarField::evaluate();
}


void ManhartWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    ManhartWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
