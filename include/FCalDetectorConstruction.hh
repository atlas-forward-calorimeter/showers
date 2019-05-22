/*
This is the header file of detctor constuction

Author: Zhaoyuan Cui
        2016 Summer
        Physics department, The University of Arizona

-

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#ifndef FCalDetectorConstruction_hh
#define FCalDetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

#include <vector>


class G4VisAttributes;
class G4LogicalVolume;


class FCalDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        // Constructor and destructor.
        FCalDetectorConstruction();
        virtual ~FCalDetectorConstruction();

        // Make everything and return the world volume.
        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();

    private:
        G4LogicalVolume *fMatrix;
        G4LogicalVolume *gapLogical;
        G4bool checkOverlaps;
        G4RotationMatrix *hexRot;

        void ConstructMaterials();  // Create materials and elements.
        void SetupGeometry();  // Create all of the geometry.

        // World logical and physical volumes.
        G4LogicalVolume*   fpWorldLogical;
        G4VPhysicalVolume* fpWorldPhysical;

        G4double height;
        G4double xStd;
};

#endif
