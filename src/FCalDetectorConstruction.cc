/*
This is the source code of the detector constuction.

Author: Zhaoyuan Cui
		2016 Summer
		Physics department, The University of Arizona

This code generates the geometry of FCal1, including liguid argon gaps and rods.
Arrangement of the rods is accomplished by following the actual region 
arrangement. Sensitive detector should be provided by the user for his or her 
own purpose.

-

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#include "FCalDetectorConstruction.hh"
#include "FCalEmCalorimeterSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4UserLimits.hh"

#include <map>


//// Constructor and Destructor

FCalDetectorConstruction::FCalDetectorConstruction() 
	: fpWorldLogical(0), fpWorldPhysical(0) {}

FCalDetectorConstruction::~FCalDetectorConstruction() {}


// Instantiate materials and volumes (once in master thread).
G4VPhysicalVolume* FCalDetectorConstruction::Construct()
{
	ConstructMaterials();
	SetupGeometry();
	return fpWorldPhysical;
}


///|////////////////////////////////////////////////////////////////////////////
//|| Sensitive Detectors
///|////////////////////////////////////////////////////////////////////////////
void FCalDetectorConstruction::ConstructSDandField()
{
	// Make a sensitive detector.
	FCalEmCalorimeterSD* aTrackerSD = new FCalEmCalorimeterSD(
		"FCalSD",               // sensitive detector name
        "FCalHitsCollection"    // detector's hits collection name
    );
	G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);

	// Make all logical volumes with the name "Gap_Logical" sensitive.
	SetSensitiveDetector("Gap_Logical", aTrackerSD, true);
}


///|////////////////////////////////////////////////////////////////////////////
//|| Materials and Elements
///|////////////////////////////////////////////////////////////////////////////
void FCalDetectorConstruction::ConstructMaterials()
{
    // For named parameters.
	G4String name, symbol;
    G4double z, Zeff;  // (effective) # of protons
    G4double a, Aeff;  // (effective) total mass
    G4double n, density, fractionmass;
	G4int ncomponents, nisotopes, natoms;

	//// Elements
    G4Element* N = new G4Element(
        name = "Nitrogen",
        symbol = "N",
        Zeff = 7.,                             
        Aeff = 14.01 * CLHEP::g / CLHEP::mole
    );
	G4Element* O = new G4Element(
		name = "Oxygen", 
        symbol = "O", 
        Zeff = 8., 
        Aeff = 16. * CLHEP::g / CLHEP::mole
	);
	G4Element* H = new G4Element(
		name = "Hydrogen", 
        symbol = "H", 
        Zeff = 1., 
        Aeff = 1.0078 * CLHEP::g / CLHEP::mole
	);
	G4Element* C = new G4Element(
		name = "Carbon", 
        symbol = "C", 
        Zeff = 6., 
        Zeff = 12 * CLHEP::g / CLHEP::mole
	);

    //// Build Materials

	//// Air
	G4Material* air = new G4Material(
		"Air", density = 1.290 * CLHEP::mg / CLHEP::cm3, ncomponents = 2
	);
	air->AddElement(N, fractionmass = 0.7);
	air->AddElement(O, fractionmass = 0.3);

	//// Vacuum
	G4Material* vacuum = new G4Material(
		"Vacuum", density = 1.e-5 * CLHEP::g / CLHEP::cm3,
		ncomponents = 1, kStateGas, CLHEP::STP_Temperature,
		2.e-2 * CLHEP::bar);
	vacuum->AddMaterial(air, fractionmass = 1.);

	//// PEEK
	G4Material* PEEK = new G4Material(
		"PEEK", density = 1.32 * CLHEP::g / CLHEP::cm3, ncomponents = 3
	);
	PEEK->AddElement(C, natoms = 19);
	PEEK->AddElement(H, natoms = 12);
	PEEK->AddElement(O, natoms = 3);

    //// Create More Materials
    new G4Material(
        name = "Copper", 
        z = 29., 
        a = 63.546 * CLHEP::g / CLHEP::mole,
        density = 8.96 * CLHEP::g / CLHEP::cm3
    );
    new G4Material(
        name = "LAr", 
        z = 18., 
        a = 39.948 * CLHEP::g / CLHEP::mole,
        density = 1.396 * CLHEP::g / CLHEP::cm3
    );
    new G4Material(
        name = "Tungsten", 
        z = 74., 
        a = 183.84 * CLHEP::g / CLHEP::mole,
        density = 19.25 * CLHEP::g / CLHEP::cm3
    );
    new G4Material(
        name = "Titanium", 
        z = 22., 
        a = 47.867 * CLHEP::g / CLHEP::mole,
        density = 4.506 * CLHEP::g / CLHEP::cm3
    );
}


///|////////////////////////////////////////////////////////////////////////////
//|| Create Detector Geometry
///|////////////////////////////////////////////////////////////////////////////
void FCalDetectorConstruction::SetupGeometry()
{
    ///|////////////////////////////////////////////////////////////////////
	//|| Dimensions and Parameters - Can be changed.
    ///|////////////////////////////////////////////////////////////////////

	//// World
	G4double worldXYZ = 200 * CLHEP::mm;

	//// Rod
	G4double rodInnerRadius = 4.5 / 2 * CLHEP::mm;
	G4double rodOuterRadius = 4.712 / 2 * CLHEP::mm;
	G4double rodMiddleZ = 10.5 / 2 * CLHEP::mm;
	G4double rodLeftRightZ = 12.25 / 2 * CLHEP::mm;
	// Total lengths of rod and tube are 35 mm.

	//// Tube
	G4double tubeInnerRadius = 5.25 / 2 * CLHEP::mm;
	G4double tubeOuterRadius = 6.26 / 2 * CLHEP::mm;
	G4double tubeMiddleZ = 8 / 2 * CLHEP::mm;
	G4double tubeLeftRightZ = 13.5 / 2 * CLHEP::mm;
	G4double tubeGap = 0.025 * CLHEP::mm;  // gap between tube parts

	//// Tube Placement
    const int numCals = 4;  // # of complete calorimeter tubes
    double tubtriside = 7.5 * CLHEP::mm;
    double tubHoleRadius = 5.8 / 2 * CLHEP::mm;

	//// Shaft
	G4double shaftRadius = 3 / 2 * CLHEP::mm;

	//// Foil and its space
	G4double foilZ = 10 / 2 * CLHEP::mm;
	G4double foilThickness = 0.02 * CLHEP::mm;
	G4double spaceThickness = 0.05 * CLHEP::mm;  // between foil & rod

	//// PEEK Box
	G4double boxX = 20. * CLHEP::mm;
	G4double boxY = 20. * CLHEP::mm;
	G4double boxZ = 43.5 / 2 * CLHEP::mm;

	//// Tungsten Plates
	G4double plateX = 30 / 2 * CLHEP::mm;
	G4double plateY = plateX;
    G4double plateZ = 3.5 / 2 * CLHEP::mm;
	G4double tunghzTot = plateX;
	G4double larGThz = 4 / 2 * CLHEP::mm;
	int tungPN = 8;		// # of front plates
	int tungBPN = 4;	// # of back plates

    //// Cryostat
    G4double cryoIIRadius = 107.95 * CLHEP::mm;       // inner wall inner radius
    G4double cryoInnerThickness = 1.91 * CLHEP::mm;   // inner wall thickness
    G4double cryoOORadius = 125 * CLHEP::mm;          // outer wall outer radius
    G4double cryoOuterThickness = 2.29 * CLHEP::mm;   // outer wall thickness
    G4double cryoPosX = -38.5 * CLHEP::mm;
    G4double cryoPosY = 0;
    // separation between inner cryostat and front plate
    G4double cryoFrontSepZ = 24.27 * CLHEP::mm;
    G4double cryoHalfLength = 100 / 2 * CLHEP::mm;
    G4double cryoStartAngle = 30 * CLHEP::degree;     // angular bounds
    G4double cryoStopAngle = 100 * CLHEP::degree;

    ///| End of Dimensions and Parameters //////////////////////////////////

    // Get materials defined in `ConstructMaterials`.
    G4Material* air = G4Material::GetMaterial("Air");
    G4Material* vacuum = G4Material::GetMaterial("Vacuum");
    G4Material* copper = G4Material::GetMaterial("Copper");
    G4Material* titanium = G4Material::GetMaterial("Titanium");
    G4Material* lar = G4Material::GetMaterial("LAr");
    G4Material* tungsten = G4Material::GetMaterial("Tungsten");
    G4Material* PEEK = G4Material::GetMaterial("PEEK");

    // Get materials from the NIST database.
    G4Material* stainlessSteel = \
        G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    //// World
	G4Box* worldSolid = new G4Box("World_Solid", worldXYZ, worldXYZ, worldXYZ);
	fpWorldLogical = new G4LogicalVolume(worldSolid, air, "World_Logical");
	fpWorldPhysical = new G4PVPlacement(
		0,					// Rotation matrix pointer
		G4ThreeVector(),	// Translation vector
		fpWorldLogical,		// Logical volume
		"World_Physical",	// Name
		0,					// Mother volume
		false,				// Unused boolean parameter
		0					// Copy number
	);

    ///|////////////////////////////////////////////////////////////////////
    //|| PEEK Box
    ///|////////////////////////////////////////////////////////////////////
    double tshifty[numCals];
    double tshiftx[numCals];
    double tshiftz[numCals];    // Zero.
    G4RotationMatrix Rth;       // No rotation.
    G4Transform3D Trh[numCals];

    for (int it = 0; it < numCals; it++)
    {
        tshiftx[it] = tubtriside / 4 * (it * 2 - 3);
        tshifty[it] = -tubtriside / 4 * sqrt(3.) * pow(-1, it);
        tshiftz[it] = 0;
        Trh[it] = G4Transform3D(
            Rth, G4ThreeVector(tshiftx[it], tshifty[it], tshiftz[it]));
    }

    G4Box* outBox = new G4Box(
        "outBox",			// Name
        boxX,				// half x thick
        boxY,				// half y thick
        boxZ				// half z thick
    );
    G4Tubs* tubHole = new G4Tubs(
        "tubHole",			// Name
        0,					// Inner radius
        tubHoleRadius,		// Outer radius
        boxZ,				// Half length in z
        0,					// Starting phi angle
        360					// Segment angle	
    );

    G4SubtractionSolid *peekBox0 = new G4SubtractionSolid(
        "peekBox0", outBox, tubHole, Trh[0]);
    G4SubtractionSolid *peekBox1 = new G4SubtractionSolid(
        "peekBox1", peekBox0, tubHole, Trh[1]);
    G4SubtractionSolid *peekBox2 = new G4SubtractionSolid(
        "peekBox2", peekBox1, tubHole, Trh[2]);
    G4SubtractionSolid *peekBox3 = new G4SubtractionSolid(
        "peekBox3", peekBox2, tubHole, Trh[3]);

    G4LogicalVolume* peekBoxLogical = new G4LogicalVolume(
        peekBox3,	// Solid
        PEEK,		// Material
        "peekBox_Logical"
    );

    new G4PVPlacement(
        0,                          // Rotation matrix
        G4ThreeVector(0, 0, 0),		// Translation vector
        peekBoxLogical,				// Logical volume
        "peekBox_Physical",
        fpWorldLogical,
        false,
        0
    );

    ///|////////////////////////////////////////////////////////////////////
    //|| Calorimeter Tubes
    ///|////////////////////////////////////////////////////////////////////

	//|| Rod (The cylinder on the inside of the LAr) ///////////////////////

	//// Rod / "Wall" (Hollow)
	G4Tubs* wallSolid = new G4Tubs(
		"Wall",							 // Name
		rodInnerRadius,					 // Inner radius
		rodOuterRadius,					 // Outer radius
		2 * rodLeftRightZ + rodMiddleZ,	 // Half length in z
		0,								 // Starting phi angle
		360								 // Segment angle
	);
	G4LogicalVolume* wallLogical = new G4LogicalVolume(
		wallSolid, titanium, "Wall_Logical");

	//// Left Rod Fill
	G4Tubs* rodLeftSolid = new G4Tubs(
		"RodLeft",			// Name
		0,					// Inner radius
		rodInnerRadius,		// Outer radius
		rodLeftRightZ,      // Half length in z
		0,					// Starting phi angle
		360					// Segment angle
	);
	G4LogicalVolume* rodLeftLogical = new G4LogicalVolume(
		rodLeftSolid, titanium,	"RodLeft_Logical");

	//// Right Rod Fill
	G4Tubs* rodRightSolid = new G4Tubs(
		"RodRight",			// Name
		0,					// Inner radius
		rodInnerRadius,		// Outer radius
		rodLeftRightZ,		// Half length in z
		0,					// Starting phi angle
		360					// Segment angle
	);
	G4LogicalVolume* rodRightLogical = new G4LogicalVolume(
		rodRightSolid, titanium, "RodRight_Logical");

	//|| Shaft (Thinnest radius. For strength.) ////////////////////////////
	G4Tubs* shaftSolid = new G4Tubs(
		"Shaft",					// Name
		0,							// Inner radius
		shaftRadius,	            // Outer radius
		rodMiddleZ,		            // Half length in z
		0,							// Starting phi angle
		360							// Segment angle
	);
	G4LogicalVolume* shaftLogical = new G4LogicalVolume(
		shaftSolid, titanium, "Shaft_Logical");

	//|| Cavity (entire cavity between outer shaft and inner rod) //////////
	//|| and the foil inside (on which the Sr90 is placed) /////////////////

	G4Tubs* cavitySolidOut = new G4Tubs(
		"CavityOut",		// Name
		shaftRadius,		// Inner radius
		rodInnerRadius,     // Outer radius
		rodMiddleZ,			// Half length in z
		0,					// Starting phi angle
		360					// Segment angle
	);

	G4Tubs* foilSolid = new G4Tubs(
		"Foil",											   // Name
		rodInnerRadius - spaceThickness - foilThickness,   // Inner radius 
		rodInnerRadius - spaceThickness,				   // Outer radius
		foilZ,											   // Half length in z
		0,												   // Starting phi angle
		360												   // Segment angle
	);

	G4SubtractionSolid* cavitySolid = new G4SubtractionSolid(
		"Cavity", cavitySolidOut, foilSolid);

	G4LogicalVolume* cavityLogical = new G4LogicalVolume(
		cavitySolid, vacuum, "Cavity_Logical");
	G4LogicalVolume* foilLogical = new G4LogicalVolume(
		foilSolid, copper, "Foil_Logical");

	//|| Gap (the cylindrical LAr gap) /////////////////////////////////////
	G4Tubs* gapSolid = new G4Tubs(
		"Gap",								// Name
		rodOuterRadius,						// Inner radius
		tubeInnerRadius,					// Outer radius
		2 * rodLeftRightZ + rodMiddleZ,		// Half length in z
		0,									// Starting phi angle
		360									// Segment angle
	);
	G4LogicalVolume* gapLogical = new G4LogicalVolume(
		gapSolid, lar, "Gap_Logical");

	//|| Tube (the FCal outer tube) ////////////////////////////////////////

    //// Middle Tube
	G4Tubs* tubeSolid = new G4Tubs(
		"Tube",				// Name
		tubeInnerRadius,    // Inner radius
		tubeOuterRadius,    // Outer radius
		tubeMiddleZ,		// Half length in z
		0,					// Starting phi angle
		360					// Segment angle
	);
	G4LogicalVolume* tubeLogical = new G4LogicalVolume(
		tubeSolid, titanium, "Tube_Logical");

	//// Left Tube
	G4Tubs* tubeSolidL = new G4Tubs(
		"TubeL",			// Name
		tubeInnerRadius,	// Inner radius
		tubeOuterRadius,	// Outer radius
		tubeLeftRightZ,		// Half length in z
		0,					// Starting phi angle
		360					// Segment angle
	);
	G4LogicalVolume* tubeLLogical =	new G4LogicalVolume(
		tubeSolidL, titanium, "TubeL_Logical");

	//// Right Tube
	G4Tubs* tubeSolidR = new G4Tubs(
		"TubeR",			// Name
		tubeInnerRadius,	// Inner radius
		tubeOuterRadius,	// Outer radius  changed from 0.2878 to 0.2875
		tubeLeftRightZ,		// Half length in z changed from 1.0 to 4.0
		0,					// Starting phi angle
		360					// Segment angle
	);
	G4LogicalVolume* tubeRLogical = new G4LogicalVolume(
		tubeSolidR, titanium, "TubeR_Logical");

    //|| Assemble and Place Calorimeters ///////////////////////////////////

    G4AssemblyVolume* assemblyTube = new G4AssemblyVolume();
    G4Transform3D Tr0;
    G4RotationMatrix Ro;

    Tr0 = G4Transform3D(Ro, G4ThreeVector(0, 0, 0));
    assemblyTube->AddPlacedVolume(
        shaftLogical, Tr0);

    Tr0 = G4Transform3D(Ro, G4ThreeVector(0, 0, -(rodMiddleZ + rodLeftRightZ)));
    assemblyTube->AddPlacedVolume(
        rodLeftLogical, Tr0);

    Tr0 = G4Transform3D(Ro, G4ThreeVector(0, 0, rodMiddleZ + rodLeftRightZ));
    assemblyTube->AddPlacedVolume(
        rodRightLogical, Tr0);

    Tr0 = G4Transform3D(Ro, G4ThreeVector(0, 0, 0));
    assemblyTube->AddPlacedVolume(
        cavityLogical, Tr0);
    assemblyTube->AddPlacedVolume(
        foilLogical, Tr0);
    assemblyTube->AddPlacedVolume(
        wallLogical, Tr0);
    assemblyTube->AddPlacedVolume(
        gapLogical, Tr0);
    assemblyTube->AddPlacedVolume(
        tubeLogical, Tr0);

    Tr0 = G4Transform3D(
        Ro, G4ThreeVector(0, 0, -(tubeMiddleZ + tubeLeftRightZ + tubeGap)));
    assemblyTube->AddPlacedVolume(
        tubeLLogical, Tr0);

    Tr0 = G4Transform3D(
        Ro, G4ThreeVector(0, 0, tubeMiddleZ + tubeLeftRightZ + tubeGap));
    assemblyTube->AddPlacedVolume(
        tubeRLogical, Tr0);

    for (int it = 0; it < numCals; it++)
    {
        // four tube
        Tr0 = G4Transform3D(
            Ro, G4ThreeVector(tshiftx[it], tshifty[it], tshiftz[it]));
        assemblyTube->MakeImprint(fpWorldLogical, Tr0);
    }
/*
    ///|////////////////////////////////////////////////////////////////////
	//|| Tungsten Plate Series
    ///|////////////////////////////////////////////////////////////////////
	G4Box* tungPlate = new G4Box("tungPlate",	// Name
		plateX,									// half x thick
		plateY,									// half y thick
		tunghzTot / (double)tungPN);			// half z thick
	G4LogicalVolume* tungPlateLogical = new G4LogicalVolume(
		tungPlate, tungsten, "tungPlate_Logical");

	G4Box* larGapTp = new G4Box("larGapTp",		// Name
		boxX,									// half x thick
		boxY,									// half y thick
		larGThz);								// half z thick
	G4LogicalVolume* larGapTpLogical = new G4LogicalVolume(
		larGapTp, lar, "Gap_Logical");

	G4double LarGapcenter;
	for (int itunp = 0; itunp < tungPN - 1; itunp++)
	{
		//// Front
		new G4PVPlacement(
			0,                  // Rotation matrix
			//Translation vector
			G4ThreeVector(0, 0, 
				boxZ \
				+ tunghzTot / tungPN \
				+ (tunghzTot / tungPN + larGThz) * 2 * itunp
			),
			tungPlateLogical,   // Logical volume
			"tungPlate_Physical",
			fpWorldLogical,
			false,
			0
		);

		LarGapcenter = boxZ + tunghzTot / tungPN * 2 + larGThz \
			+ (tunghzTot / tungPN + larGThz) * 2 * itunp;

		new G4PVPlacement(
			0,										//Rotation matrix
			G4ThreeVector(0, 0, LarGapcenter),		//Translation vector
			larGapTpLogical,						//Logical volume
			"larGapTp_Physical",
			fpWorldLogical,
			false,
			0
		);

		//// Back
		if (itunp < tungBPN - 1)
		{
			new G4PVPlacement(
				0,  //Rotation matrix
				//Translation vector
				G4ThreeVector(
					0, 0, 
					-(
						boxZ
						+ tunghzTot / tungPN
						+ (tunghzTot / tungPN + larGThz) * 2 * itunp
					)
				),
				tungPlateLogical,  //Logical volume
				"tungPlate_Physical",
				fpWorldLogical,
				false,
				0
			);

			new G4PVPlacement(
				0,										//Rotation matrix
				G4ThreeVector(0, 0, -LarGapcenter),		//Translation vector
				larGapTpLogical,                        //Logical volume
				"larGapTp_Physical",
				fpWorldLogical,
				false,
				0
			);
		}
	}

	//// Last Piece - Front
	new G4PVPlacement(
		0,                  //Rotation matrix
		//Translation vector
		G4ThreeVector(0, 0, 
			boxZ \
			+ tunghzTot / tungPN \
			+ (tunghzTot / tungPN + larGThz) * 2 * (tungPN - 1)
		),
		tungPlateLogical,   //Logical volume
		"tungPlate_Physical",
		fpWorldLogical,
		false,
		0
	);

	//// Last Piece - Back
	new G4PVPlacement(
		0,                  //Rotation matrix
		//Translation vector
		G4ThreeVector(0, 0, 
			-(
				boxZ
				+ tunghzTot / tungPN
				+ (tunghzTot / tungPN + larGThz) * 2 * (tungBPN - 1)
			)
		),
		tungPlateLogical,   //Logical volume
		"tungPlate_Physical",
		fpWorldLogical,
		false,
		0
	);
    ///| End of Tungsten Plate Series //////////////////////////////////////
*/
    ///|////////////////////////////////////////////////////////////////////
    //|| Tungsten Plate Series
    ///|////////////////////////////////////////////////////////////////////
    G4Box* tungPlate = new G4Box("tungPlate",	// Name
        plateX,									// half x thick
        plateY,									// half y thick
        plateZ);			                    // half z thick
    G4LogicalVolume* tungPlateLogical = new G4LogicalVolume(
        tungPlate, tungsten, "tungPlate_Logical");

    G4Box* larGapTp = new G4Box("larGapTp",		// Name
        boxX,									// half x thick
        boxY,									// half y thick
        larGThz);								// half z thick
    G4LogicalVolume* larGapTpLogical = new G4LogicalVolume(
        larGapTp, lar, "Gap_Logical");

    G4double LarGapcenter;
    for (int itunp = 0; itunp < tungPN - 1; itunp++)
    {
        //// Front Plates and Gaps
        new G4PVPlacement(
            0,                  // Rotation matrix
            //Translation vector
            G4ThreeVector(0, 0,
                boxZ \
                + plateZ \
                + 2 * (plateZ + larGThz) * itunp
            ),
            tungPlateLogical,   // Logical volume
            "tungPlate_Physical",
            fpWorldLogical,
            false,
            0
        );

        LarGapcenter = boxZ + 2 * plateZ + larGThz \
            + 2 * (plateZ + larGThz) * itunp;

        new G4PVPlacement(
            0,										//Rotation matrix
            G4ThreeVector(0, 0, LarGapcenter),		//Translation vector
            larGapTpLogical,						//Logical volume
            "larGapTp_Physical",
            fpWorldLogical,
            false,
            0
        );

        //// Back Plates and Gaps
        if (itunp < tungBPN - 1)
        {
            new G4PVPlacement(
                0,  //Rotation matrix
                //Translation vector
                G4ThreeVector(
                    0, 0,
                    -(boxZ
                      + plateZ
                      + 2 * (plateZ + larGThz) * itunp)
                ),
                tungPlateLogical,  //Logical volume
                "tungPlate_Physical",
                fpWorldLogical,
                false,
                0
            );

            new G4PVPlacement(
                0,										//Rotation matrix
                G4ThreeVector(0, 0, -LarGapcenter),		//Translation vector
                larGapTpLogical,                        //Logical volume
                "larGapTp_Physical",
                fpWorldLogical,
                false,
                0
            );
        }
    }

    //// Last Plate (No LAr Gap) - Front
    new G4PVPlacement(
        0,                  //Rotation matrix
        //Translation vector
        G4ThreeVector(0, 0,
            boxZ \
            + plateZ \
            + 2 * (plateZ + larGThz) * (tungPN - 1)
        ),
        tungPlateLogical,   //Logical volume
        "tungPlate_Physical",
        fpWorldLogical,
        false,
        0
    );

    //// Last Plate (No LAr Gap) - Back
    new G4PVPlacement(
        0,                  //Rotation matrix
        //Translation vector
        G4ThreeVector(0, 0,
            -(boxZ
              + plateZ
              + 2 * (plateZ + larGThz) * (tungBPN - 1))
        ),
        tungPlateLogical,   //Logical volume
        "tungPlate_Physical",
        fpWorldLogical,
        false,
        0
    );
    ///| End of Tungsten Plate Series //////////////////////////////////////

    ///|////////////////////////////////////////////////////////////////////
    //|| Cryostat
    ///|////////////////////////////////////////////////////////////////////

    // Position and rotation.
    G4RotationMatrix cryoRotation;
    cryoRotation.rotateX(90 * CLHEP::degree);
    G4double cryoPosZ = \
        cryoFrontSepZ
        + boxZ + 2 * (plateZ + larGThz) * (tungPN - 1) + 2 * plateZ
        - sqrt(pow(cryoIIRadius, 2) - pow(cryoPosX, 2));
    G4Transform3D cryoTransform = G4Transform3D(
        cryoRotation,
        G4ThreeVector(cryoPosX, cryoPosY, cryoPosZ)
    );

    //// Cryostat Inner Wall
    G4Tubs* cryoInnerWallSolid = new G4Tubs(
        "Cryo_Inner_Wall",                  // Name
        cryoIIRadius,	                    // Inner radius
        cryoIIRadius + cryoInnerThickness,	// Outer radius
        cryoHalfLength,                     // Half length
        cryoStartAngle,				        // Starting phi angle
        cryoStopAngle - cryoStartAngle      // Segment angle
    );
    G4LogicalVolume* cryoInnerWallLogical = new G4LogicalVolume(
        cryoInnerWallSolid, stainlessSteel, "Cryo_Inner_Wall_Logical");
    new G4PVPlacement(
        cryoTransform,
        cryoInnerWallLogical,
        "Cryo_Inner_Wall_Physical",
        fpWorldLogical,
        false,
        0
    );

    //// Cryostat Middle, Liquid Argon
    G4Tubs* cryoMiddleSolid = new G4Tubs(
        "Cryo_Middle",                      // Name
        cryoIIRadius + cryoInnerThickness,	// Inner radius
        cryoOORadius - cryoOuterThickness,	// Outer radius
        cryoHalfLength,                     // Half length
        cryoStartAngle,					    // Starting phi angle
        cryoStopAngle - cryoStartAngle		// Segment angle
    );
    G4LogicalVolume* cryoMiddleLogical = new G4LogicalVolume(
        cryoMiddleSolid, lar, "Cryo_Middle_Logical");
    new G4PVPlacement(
        cryoTransform,
        cryoMiddleLogical,
        "Cryo_Middle_Physical",
        fpWorldLogical,
        false,
        0
    );

    //// Cryostat Outer Wall
    G4Tubs* cryoOuterWallSolid = new G4Tubs(
        "Cryo_Inner_Wall",                  // Name
        cryoOORadius - cryoOuterThickness,	// Inner radius
        cryoOORadius,	                    // Outer radius
        cryoHalfLength,                     // Half length
        cryoStartAngle,					    // Starting phi angle
        cryoStopAngle - cryoStartAngle		// Segment angle
    );
    G4LogicalVolume* cryoOuterWallLogical = new G4LogicalVolume(
        cryoOuterWallSolid, stainlessSteel, "Cryo_Outer_Wall_Logical");
    new G4PVPlacement(
        cryoTransform,
        cryoOuterWallLogical,
        "Cryo_Outer_Wall_Physical",
        fpWorldLogical,
        false,
        0
    );

	//|| Visualization /////////////////////////////////////////////////////

	//// Change color palette here
	G4VisAttributes* colors[] = {
		new G4VisAttributes(G4Colour(1.0, 0.2, 0.2, 1.0)),	// Hot Red      1.0
		new G4VisAttributes(G4Colour(0.0, 0.5, 0.5, 1.0)),	// Argon Green  0.15
		new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 1.0)),	// Silver       0.3
		new G4VisAttributes(G4Colour(0.4, 0.15, 0.4, 1.0))	// Tungsten     0.3
	};
	for (const auto& color : colors) {
		color->SetForceSolid(true);
	}

	//// Do art here
	std::map<G4LogicalVolume*, G4VisAttributes*> colorMap = {
		{foilLogical, colors[0]},
		{gapLogical, colors[1]},
		{larGapTpLogical, colors[1]},
		{tungPlateLogical, colors[3]},
		{shaftLogical, colors[2]},
		{wallLogical, colors[2]},
		{tubeLogical, colors[2]},
        {cryoInnerWallLogical, colors[2]},
        {cryoMiddleLogical, colors[1]},
        {cryoOuterWallLogical, colors[2]},
	};
	for (const auto& p : colorMap) {
		p.first->SetVisAttributes(p.second);
	}
	///| End Of Visualization //////////////////////////////////////////////

	// Tracking length in each volume:

	G4double maxStep = 0.001*CLHEP::mm; // 1 mu tracking
	G4UserLimits * limits = new G4UserLimits(maxStep);

	fpWorldLogical->SetUserLimits(limits);
	shaftLogical->SetUserLimits(limits);
	cavityLogical->SetUserLimits(limits);
	foilLogical->SetUserLimits(limits);
	wallLogical->SetUserLimits(limits);
	gapLogical->SetUserLimits(limits);
	tubeLogical->SetUserLimits(limits);
}