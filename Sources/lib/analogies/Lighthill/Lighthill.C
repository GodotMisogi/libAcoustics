/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "Lighthill.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

#include "fvc.H"
#include "sampledPatch.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Lighthill, 0);

    addToRunTimeSelectionTable(functionObject, Lighthill, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Lighthill::Lighthill
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    AcousticAnalogy
    (
        name,
        runTime,
        dict
    ),
    c_(vector::zero),
    F_(vector::zero, obr_.time().value())
{
    this->read(dict);
    F_.resize(1);
}

Foam::functionObjects::Lighthill::Lighthill
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    AcousticAnalogy
    (
        name,
        obr,
        dict
    ),
    c_(vector::zero),
    F_(vector::zero, obr_.time().value())
{
    this->read(dict);
    F_.resize(1);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Lighthill::~Lighthill()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Lighthill::read(const dictionary& dict)
{
    if (!AcousticAnalogy::read(dict))
    {
        return false;
    }

    calcDistances();

    return true;
}

void Foam::functionObjects::Lighthill::calcDistances()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    vectorField ci(0);
    scalar ni = 0;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        ci.append(mesh.boundary()[patchi].Cf());
        ni += scalar(ci.size());
    }
    reduce (ni, sumOp<scalar>());
    c_ = gSum(ci) / ni;
}

bool Foam::functionObjects::Lighthill::execute()
{
    return true;
}

bool Foam::functionObjects::Lighthill::write()
{
    //use forces library to calculate forces acting on patch
    if (!AcousticAnalogy::write())
    {
        return false;
    }

    return true;
}

void Foam::functionObjects::Lighthill::correct()
{
    calcForcesMoment();

    F_.value(0) = forceEff();
    vector dotF = F_.dot(obr_.time().value(), 0);

    if (Pstream::master() || !Pstream::parRun())
    {
        scalar coeff1_3d = 1. / 4. / Foam::constant::mathematical::pi / c0_;

//        scalar coeff1_2d = 1. / 2.828427 / Foam::constant::mathematical::pi / sqrt(c0_);
//        scalar coeff2_2d = 1. / 2. / Foam::constant::mathematical::pi;
//        scalar coeff3_2d = sqrt(c0_) / 5.656854 / Foam::constant::mathematical::pi;

        forAll (observers_, iObs)
        {
            SoundObserver& obs = observers_[iObs];
            vector l = obs.position() - c_;
            scalar r = mag(l);
            scalar oap = l & (dotF + c0_ * F_.value(0) / r) * coeff1_3d / r / r;;
            if (dRef_ > 0.0)
            {
                oap /= dRef_;
            }
            obs.apressure(oap);
        }
    }
}


// ************************************************************************* //


Foam::tmp<Foam::vectorField> Foam::functionObjects::Lighthill::surfaceStressDivergence(const sampledSurface& surface) const
{
    tmp<Field<vector>> divTSampled;
    
    volScalarField p (obr_.lookupObject<volScalarField>("p"));
    volVectorField U (obr_.lookupObject<volVectorField>("U"));

    volVectorField gradp = fvc::grad(p);

    volTensorField momConv = this->rho()() * U * U;
    volVectorField divMomConv = fvc::div(momConv);

    divTSampled = sampleOrInterpolate<vector>(gradp,surface);
    if (p.dimensions() != dimPressure)
    {
        divTSampled.ref() *= rhoRef_;
    }

    divTSampled.ref() += sampleOrInterpolate<vector>(divMomConv,surface);

    return divTSampled;
}
