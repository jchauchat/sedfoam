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

#include "Camenen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Camenen, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        Camenen,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Camenen::Camenen
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    dragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Camenen::~Camenen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Camenen::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta(scalar(1) - alpha_);

// find a way to get access to gravity vector g
    dimensionedScalar magg
    (
        "magg",
        dimensionSet(0, 1, -2, 0, 0, 0, 0),
        scalar(9.81e0)
    );
    volScalarField phi((phasea_.rho()-phaseb_.rho())/
                        (phasea_.rhoFloc()-phaseb_.rho())*alpha_
                      );
    volScalarField arg1(phasea_.Ws0()*
                        pow(beta, phasea_.dimFrac()/2)*
                        pow(max(1-phi, 1e-6), phasea_.dimFrac()/2-1)*
                        pow(max(1-phi/phasea_.phiMax(), 1e-6), phasea_.phiMax())
                       );
    volScalarField arg2(phasea_.WsGel()*
                        pow(max(phasea_.Xi()*alpha_/phasea_.alphaGel(), 1e-6),
                        1-2/(3-phasea_.dimFrac()))
                       );

    return (phasea_.rho()-phaseb_.rho())*magg
            /( neg(alpha_-phasea_.alphaGel()/phasea_.Xi()) *beta* arg1
            +pos(alpha_-phasea_.alphaGel()/phasea_.Xi()) *beta* arg2
           );
}

// ************************************************************************* //
