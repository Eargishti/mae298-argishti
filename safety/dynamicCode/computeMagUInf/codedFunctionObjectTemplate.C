/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedFunctionObjectTemplate.H"
#include "volFields.H"
#include "read.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(computeMagUInfFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    computeMagUInfFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 3c67dff08d2577d0ec2a65ce9a683703864de19c
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void computeMagUInf_3c67dff08d2577d0ec2a65ce9a683703864de19c(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& computeMagUInfFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

computeMagUInfFunctionObject::computeMagUInfFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

computeMagUInfFunctionObject::~computeMagUInfFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool computeMagUInfFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read computeMagUInf sha1: 3c67dff08d2577d0ec2a65ce9a683703864de19c\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


Foam::wordList computeMagUInfFunctionObject::fields() const
{
    if (false)
    {
        Info<<"fields computeMagUInf sha1: 3c67dff08d2577d0ec2a65ce9a683703864de19c\n";
    }

    wordList fields;
//{{{ begin code
    
//}}} end code

    return fields;
}


bool computeMagUInfFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute computeMagUInf sha1: 3c67dff08d2577d0ec2a65ce9a683703864de19c\n";
    }

//{{{ begin code
    #line 58 "/home/edgar/OpenFOAM-dev/tutorials/incompressibleFluid/EdgarFoil/system/functions/computeMagUInf"

            const volVectorField& U = mesh().lookupObject<volVectorField>("U");

            vector Umean = gSum(U.internalField()*mesh().V()) / gSum(mesh().V());
            scalar Umag = mag(Umean);

            Info<< "Computed freestream velocity magnitude Uinf = " << Umag << endl;

            OFstream file("constant/computedMagUInf");
            file << Umag << endl;
        
//}}} end code

    return true;
}


bool computeMagUInfFunctionObject::write()
{
    if (false)
    {
        Info<<"write computeMagUInf sha1: 3c67dff08d2577d0ec2a65ce9a683703864de19c\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool computeMagUInfFunctionObject::end()
{
    if (false)
    {
        Info<<"end computeMagUInf sha1: 3c67dff08d2577d0ec2a65ce9a683703864de19c\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

