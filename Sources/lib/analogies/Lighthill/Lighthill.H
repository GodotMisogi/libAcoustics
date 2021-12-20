/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::Lighthill

Group
    grpLighthillFunctionObjects

Description


SeeAlso


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef Lighthill_H
#define Lighthill_H

#include "AcousticAnalogy.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class Lighthill Declaration
\*---------------------------------------------------------------------------*/

class Lighthill
:
    public AcousticAnalogy
{

        //- Center of surface for sampling
        vector c_;

        //- Force acting on patch
        tdPrimitive<vector> F_;

protected:

        //- Disallow default bitwise copy construct
        Lighthill(const Lighthill&);

        //- Disallow default bitwise assignment
        void operator=(const Lighthill&);

        //-
        void calcDistances();

        //-
        void correct();

        //- Interpolation scheme to use for sampling
        word interpolationScheme_;

        //- Sample and write a particular volume field
        template<class Type>
        tmp<Field<Type> > sampleOrInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const sampledSurface& surface
        ) const;

        tmp<vectorField> surfaceStressDivergence(const sampledSurface& surface) const;

public:

    //- Runtime type information
    TypeName("Lighthill");


    // Constructors

        //- Construct from given Time and dictionary.
        Lighthill
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        Lighthill
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~Lighthill();


    // Member Functions

        //- Read the Lighthill data
        virtual bool read(const dictionary&);

        //-
        virtual bool execute();

        //- Write the Lighthill
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LighthillTemplates.C"

#endif
// ************************************************************************* //
