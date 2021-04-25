/*
Global definitions.

*/

#ifndef FCalGlobals_h
#define FCalGlobals_h 1

namespace FCal {
    const G4double rodInnerRadius = 4.5 / 2 * CLHEP::mm;
    const G4double spaceThickness = 0.05 * CLHEP::mm;  // between foil & rod (0.05 original)
    const G4double foilOuterRadius = rodInnerRadius - spaceThickness;
    const G4double foilZ = 10 / 2 * CLHEP::mm;  // half width
    //const G4double rodOuterRadius = 4.712 / 2 * CLHEP::mm;  // Rod
    //const G4double tubeInnerRadius = 5.25 / 2 * CLHEP::mm;  // Tube
}

#endif
