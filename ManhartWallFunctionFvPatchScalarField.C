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
    const fvPatchVectorField Gbpr = gradp.boundaryField()[patchi];	
   // const  fvPatchScalarField& bpr = pr.boundaryField()[patchi]; 
    
   //  volVectorField Gbpr = fvc::grad(bpr);
    scalarField GradP = Gbpr.component(0); 

    forAll(nuSgsw, facei) // facei for how many faces ? = boundary field?
    {
        scalar magUpara = magUp[facei];  // Wall Parallel "U" velocity at the cell adjacent to the wall

        scalar Pgrad = GradP[facei];
       
        std::cout << "\n";
        std::cout << "magUpara is \t" << magUpara ;       
       
        scalar tau = (nuSgsw[facei] + nuw[facei])*magFaceGradU[facei];  // Wall shear stress (tau) definition.        
       // scalar utau = sqrt((nuSgsw[facei] + nuw[facei])*magFaceGradU[facei]);
        std::cout << "magFaceGradU[facei] is \t" << magFaceGradU[facei]; 
        
        scalar utau = sqrt(mag(tau)); // Friction Velocity (Utau) definition. 

        scalar uP = pow(mag(( nuSgsw[facei] +   nuw[facei])*Pgrad),(1.0/3.0));  // Defining streamwise pressure based velocity (nuw or nusgs?)

        std::cout <<  "GradP[facei] is \t" << Pgrad << "\n"; //DEBUG
        std::cout << "nuSgsw is \t" << nuSgsw[facei] << "\t and nuw is \t" << nuw[facei]; // DEBUG         

        scalar utauP = sqrt(sqr(utau) + sqr(uP)); 
 
        std::cout << "utauP is \t" << utauP << "\n"; //DEBUG
        std::cout  << "utau is" << utau << "\t and uP is \t" << uP << "\n" ;  // DEBUG 

        scalar alpha = sqr(utau)/sqr(utauP);

        std::cout <<  "Alpha is \t" << alpha << "\n";

     //   scalar Ustar = magUpara/utauP;

  //      scalar  ystar =  utauP/(ry[facei]*nuw[facei]);   
        std::cout << " Newton Raphson Loop Start " << "\n";
        
        if (utau > VSMALL)   
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {  //Newton-Raphson Solution
               

               scalar YoverNu =  1.0/(ry[facei]*nuw[facei]);
               scalar utauPsqr = sqr(utau) + sqr(uP);
               // F
               scalar f =  sign(tau)*(sqr(utau)/sqrt(utauPsqr))*YoverNu 
                           + sign(Pgrad)*0.5*(pow(uP,3)/sqrt(utauPsqr))*sqr(YoverNu) 
                           - magUpara/sqrt(utauPsqr);
               
               // DEBUG
               std::cout << " The new uP is" << uP << "\n";
               std::cout << "The prev. iteration  utauP is" << sqrt(utauPsqr) << "\t The prev. iter utau is" << utau << "\n"; //DEBUG
               std::cout << "YoverNu is \t" << YoverNu << "\n"; //DEBUG
               // F DERIVATIVE
               scalar df =  sign(tau)*(2*utau/sqrt(utauPsqr) - pow(utau,3)/pow(utauPsqr,1.5))*YoverNu 
                            + magUpara*(utau/pow(utauPsqr,1.5))
                            - sign(Pgrad)*0.5*(pow(uP,3)*(1.0/pow(utauPsqr,1.5))*utau)*sqr(YoverNu);

               //DEBUG
               std::cout << " NR Values " << "\n";     //DEBUG
               std::cout << " the 'f' is " << f << "\n";
               std::cout << " the 'df' is " << df <<  "\n";
      
              
               // correction   
               scalar utauNew = utau - f/df;
 
               err = mag((utau - utauNew)/utau);

               utau = utauNew;
               std::cout << "utauNew is" << utau << "\n";
               std::cout << "Error is" << err << "\n";

               tau = sqr(utau); // sqr(max(utau,0)); 
             
            } while (utau > VSMALL && err > 0.01 && ++iter < 10);        // Setting tolerance criteria for Utau and max limit of NR iterations
          

             std::cout << " New tau is \t" << tau; //DEBUG
             std::cout << " The new magFaceGradU is \t" << magFaceGradU[facei];
             std::cout << " The new nuSgs is \t" << nuSgsw[facei] << " \t and the new nuw[facei] \t" << nuw[facei];

             scalar nuCorr = (tau/magFaceGradU[facei])*(1 + (nuSgsw[facei]/nuw[facei])) - nuw[facei] - nuSgsw[facei];
             
             std::cout << "nuCorr is \t" << nuCorr;  //DEBUG

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
