// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS contributors $
//! Element and isotope data transcribed from OpenMS4-core ElementDB.cpp at 7c029e8.
//! Atomic masses and abundances preserve the declared upstream tables.
//! Iridium deliberately uses its declared iridium table, correcting the upstream
//! buildElement_ call that accidentally passes the rhenium table.

// Keep the upstream decimal literals intact for source-level auditability.
#![allow(clippy::excessive_precision)]

use super::{Element, Isotope};

#[rustfmt::skip]
pub(super) static ELEMENTS: &[Element] = &[
    Element { name: "Hydrogen", symbol: "H", atomic_number: 1, isotopes: &[
        Isotope { mass_number: 1, mass: 1.0078250319, abundance: 0.999885 },
        Isotope { mass_number: 2, mass: 2.01410178, abundance: 0.000115 },
        Isotope { mass_number: 3, mass: 3.01604927, abundance: 0.0 },
    ] },
    Element { name: "Helium", symbol: "He", atomic_number: 2, isotopes: &[
        Isotope { mass_number: 3, mass: 3.0160293191, abundance: 1.34e-06 },
        Isotope { mass_number: 4, mass: 4.00260325415, abundance: 0.9999986599999999 },
    ] },
    Element { name: "Lithium", symbol: "Li", atomic_number: 3, isotopes: &[
        Isotope { mass_number: 6, mass: 6.015122, abundance: 0.0759 },
        Isotope { mass_number: 7, mass: 7.016004, abundance: 0.9240999999999999 },
    ] },
    Element { name: "Beryllium", symbol: "Be", atomic_number: 4, isotopes: &[
        Isotope { mass_number: 9, mass: 9.0121822, abundance: 1.0 },
    ] },
    Element { name: "Boron", symbol: "B", atomic_number: 5, isotopes: &[
        Isotope { mass_number: 10, mass: 10.012937000000001, abundance: 0.19899999999999998 },
        Isotope { mass_number: 11, mass: 11.009304999999999, abundance: 0.8009999999999999 },
    ] },
    Element { name: "Carbon", symbol: "C", atomic_number: 6, isotopes: &[
        Isotope { mass_number: 12, mass: 12.0, abundance: 0.9893000000000001 },
        Isotope { mass_number: 13, mass: 13.003355000000001, abundance: 0.010700000000000001 },
    ] },
    Element { name: "Nitrogen", symbol: "N", atomic_number: 7, isotopes: &[
        Isotope { mass_number: 14, mass: 14.003074, abundance: 0.9963200000000001 },
        Isotope { mass_number: 15, mass: 15.000109, abundance: 0.00368 },
    ] },
    Element { name: "Oxygen", symbol: "O", atomic_number: 8, isotopes: &[
        Isotope { mass_number: 16, mass: 15.994915000000001, abundance: 0.9975700000000001 },
        Isotope { mass_number: 17, mass: 16.999132, abundance: 0.00037999999999999997 },
        Isotope { mass_number: 18, mass: 17.999168999999998, abundance: 0.0020499999999999997 },
    ] },
    Element { name: "Fluorine", symbol: "F", atomic_number: 9, isotopes: &[
        Isotope { mass_number: 19, mass: 18.99840322, abundance: 1.0 },
    ] },
    Element { name: "Neon", symbol: "Ne", atomic_number: 10, isotopes: &[
        Isotope { mass_number: 20, mass: 19.99244018, abundance: 0.9048 },
        Isotope { mass_number: 21, mass: 20.9938467, abundance: 0.0027 },
        Isotope { mass_number: 22, mass: 21.9913851, abundance: 0.0925 },
    ] },
    Element { name: "Sodium", symbol: "Na", atomic_number: 11, isotopes: &[
        Isotope { mass_number: 23, mass: 22.989769280899999, abundance: 1.0 },
    ] },
    Element { name: "Magnesium", symbol: "Mg", atomic_number: 12, isotopes: &[
        Isotope { mass_number: 24, mass: 23.985042, abundance: 0.7898999999999999 },
        Isotope { mass_number: 25, mass: 24.985837, abundance: 0.1 },
        Isotope { mass_number: 26, mass: 25.982593000000001, abundance: 0.1101 },
    ] },
    Element { name: "Aluminium", symbol: "Al", atomic_number: 13, isotopes: &[
        Isotope { mass_number: 27, mass: 26.981538629999999, abundance: 1.0 },
    ] },
    Element { name: "Silicon", symbol: "Si", atomic_number: 14, isotopes: &[
        Isotope { mass_number: 28, mass: 27.976926532499999, abundance: 0.9220999999999999 },
        Isotope { mass_number: 29, mass: 28.9764947, abundance: 0.0467 },
        Isotope { mass_number: 30, mass: 29.973770170000002, abundance: 0.031 },
    ] },
    Element { name: "Phosphorus", symbol: "P", atomic_number: 15, isotopes: &[
        Isotope { mass_number: 31, mass: 30.973761490000001, abundance: 1.0 },
    ] },
    Element { name: "Sulfur", symbol: "S", atomic_number: 16, isotopes: &[
        Isotope { mass_number: 32, mass: 31.972070729999999, abundance: 0.9493 },
        Isotope { mass_number: 33, mass: 32.971457999999998, abundance: 0.0076 },
        Isotope { mass_number: 34, mass: 33.967866999999998, abundance: 0.0429 },
        Isotope { mass_number: 36, mass: 35.967081, abundance: 0.0002 },
    ] },
    Element { name: "Chlorine", symbol: "Cl", atomic_number: 17, isotopes: &[
        Isotope { mass_number: 35, mass: 34.968852679999998, abundance: 0.7576 },
        Isotope { mass_number: 37, mass: 36.965902589999999, abundance: 0.24239999999999998 },
    ] },
    Element { name: "Argon", symbol: "Ar", atomic_number: 18, isotopes: &[
        Isotope { mass_number: 36, mass: 35.967545106000003, abundance: 0.003336 },
        Isotope { mass_number: 38, mass: 37.9627324, abundance: 0.000629 },
        Isotope { mass_number: 40, mass: 39.9623831225, abundance: 0.996035 },
    ] },
    Element { name: "Potassium", symbol: "K", atomic_number: 19, isotopes: &[
        Isotope { mass_number: 39, mass: 38.963706680000001, abundance: 0.932581 },
        Isotope { mass_number: 40, mass: 39.963998480000001, abundance: 0.000117 },
        Isotope { mass_number: 41, mass: 40.961825760000004, abundance: 0.067302 },
    ] },
    Element { name: "Calcium", symbol: "Ca", atomic_number: 20, isotopes: &[
        Isotope { mass_number: 40, mass: 39.962590980000002, abundance: 0.96941 },
        Isotope { mass_number: 42, mass: 41.958618010000002, abundance: 0.00647 },
        Isotope { mass_number: 43, mass: 42.958766599999997, abundance: 0.00135 },
        Isotope { mass_number: 44, mass: 43.955481800000001, abundance: 0.02086 },
        Isotope { mass_number: 46, mass: 45.953692599999997, abundance: 4e-05 },
        Isotope { mass_number: 48, mass: 47.952534, abundance: 0.00187 },
    ] },
    Element { name: "Scandium", symbol: "Sc", atomic_number: 21, isotopes: &[
        Isotope { mass_number: 45, mass: 44.955910000000003, abundance: 1.0 },
    ] },
    Element { name: "Titanium", symbol: "Ti", atomic_number: 22, isotopes: &[
        Isotope { mass_number: 46, mass: 45.952631599999997, abundance: 0.0825 },
        Isotope { mass_number: 47, mass: 46.951763100000001, abundance: 0.07440000000000001 },
        Isotope { mass_number: 48, mass: 47.947946299999998, abundance: 0.7372 },
        Isotope { mass_number: 49, mass: 48.947870000000002, abundance: 0.0541 },
        Isotope { mass_number: 50, mass: 49.944791199999997, abundance: 0.0518 },
    ] },
    Element { name: "Vanadium", symbol: "V", atomic_number: 23, isotopes: &[
        Isotope { mass_number: 50, mass: 49.947158500000001, abundance: 0.0025 },
        Isotope { mass_number: 51, mass: 50.943959499999998, abundance: 0.9975 },
    ] },
    Element { name: "Chromium", symbol: "Cr", atomic_number: 24, isotopes: &[
        Isotope { mass_number: 50, mass: 49.946044200000003, abundance: 0.043449999999999996 },
        Isotope { mass_number: 52, mass: 51.940507500000003, abundance: 0.83789 },
        Isotope { mass_number: 53, mass: 52.940649399999998, abundance: 0.09501 },
        Isotope { mass_number: 54, mass: 53.938880400000002, abundance: 0.02365 },
    ] },
    Element { name: "Manganese", symbol: "Mn", atomic_number: 25, isotopes: &[
        Isotope { mass_number: 55, mass: 54.938049999999997, abundance: 1.0 },
    ] },
    Element { name: "Ferrum", symbol: "Fe", atomic_number: 26, isotopes: &[
        Isotope { mass_number: 54, mass: 53.939610500000001, abundance: 0.058449999999999995 },
        Isotope { mass_number: 56, mass: 55.934937499999997, abundance: 0.91754 },
        Isotope { mass_number: 57, mass: 56.935394000000002, abundance: 0.021191 },
        Isotope { mass_number: 58, mass: 57.933275600000002, abundance: 0.002819 },
    ] },
    Element { name: "Cobalt", symbol: "Co", atomic_number: 27, isotopes: &[
        Isotope { mass_number: 59, mass: 58.933194999999998, abundance: 1.0 },
    ] },
    Element { name: "Nickel", symbol: "Ni", atomic_number: 28, isotopes: &[
        Isotope { mass_number: 58, mass: 57.935347999999998, abundance: 0.680169 },
        Isotope { mass_number: 60, mass: 59.930790999999999, abundance: 0.262231 },
        Isotope { mass_number: 61, mass: 60.931060000000002, abundance: 0.011399 },
        Isotope { mass_number: 62, mass: 61.928348999999997, abundance: 0.036345 },
        Isotope { mass_number: 64, mass: 63.927970000000002, abundance: 0.009256 },
    ] },
    Element { name: "Copper", symbol: "Cu", atomic_number: 29, isotopes: &[
        Isotope { mass_number: 63, mass: 62.929600999999998, abundance: 0.6917 },
        Isotope { mass_number: 65, mass: 64.927794000000006, abundance: 0.30829999999999996 },
    ] },
    Element { name: "Zinc", symbol: "Zn", atomic_number: 30, isotopes: &[
        Isotope { mass_number: 64, mass: 63.929147, abundance: 0.4863 },
        Isotope { mass_number: 66, mass: 65.926036999999994, abundance: 0.27899999999999997 },
        Isotope { mass_number: 67, mass: 66.927131000000003, abundance: 0.040999999999999995 },
        Isotope { mass_number: 68, mass: 67.924847999999997, abundance: 0.1875 },
        Isotope { mass_number: 70, mass: 69.925325000000001, abundance: 0.0062 },
    ] },
    Element { name: "Gallium", symbol: "Ga", atomic_number: 31, isotopes: &[
        Isotope { mass_number: 69, mass: 68.925573600000007, abundance: 0.60108 },
        Isotope { mass_number: 71, mass: 70.924701299999995, abundance: 0.39892000000000005 },
    ] },
    Element { name: "Germanium", symbol: "Ge", atomic_number: 32, isotopes: &[
        Isotope { mass_number: 70, mass: 69.924247399999999, abundance: 0.20379999999999998 },
        Isotope { mass_number: 72, mass: 71.922075800000002, abundance: 0.2731 },
        Isotope { mass_number: 73, mass: 72.9234589, abundance: 0.0776 },
        Isotope { mass_number: 74, mass: 73.921177799999995, abundance: 0.36719999999999997 },
        Isotope { mass_number: 76, mass: 75.921401, abundance: 0.0776 },
    ] },
    Element { name: "Arsenic", symbol: "As", atomic_number: 33, isotopes: &[
        Isotope { mass_number: 75, mass: 74.921596500000007, abundance: 1.0 },
    ] },
    Element { name: "Selenium", symbol: "Se", atomic_number: 34, isotopes: &[
        Isotope { mass_number: 74, mass: 73.922476399999994, abundance: 0.00889 },
        Isotope { mass_number: 76, mass: 75.919213600000006, abundance: 0.09366 },
        Isotope { mass_number: 77, mass: 76.919914000000006, abundance: 0.07635 },
        Isotope { mass_number: 78, mass: 77.917309099999997, abundance: 0.23772 },
        Isotope { mass_number: 80, mass: 79.916521299999999, abundance: 0.49607 },
        Isotope { mass_number: 82, mass: 81.916699399999999, abundance: 0.08731 },
    ] },
    Element { name: "Bromine", symbol: "Br", atomic_number: 35, isotopes: &[
        Isotope { mass_number: 79, mass: 78.918337100000002, abundance: 0.5069 },
        Isotope { mass_number: 81, mass: 80.916290599999996, abundance: 0.49310000000000004 },
    ] },
    Element { name: "Krypton", symbol: "Kr", atomic_number: 36, isotopes: &[
        Isotope { mass_number: 78, mass: 77.920400000000001, abundance: 0.0034999999999999996 },
        Isotope { mass_number: 80, mass: 79.916380000000004, abundance: 0.0225 },
        Isotope { mass_number: 82, mass: 81.913482000000002, abundance: 0.11599999999999999 },
        Isotope { mass_number: 83, mass: 82.914135000000002, abundance: 0.115 },
        Isotope { mass_number: 84, mass: 83.911507, abundance: 0.57 },
        Isotope { mass_number: 86, mass: 85.910616000000005, abundance: 0.17300000000000001 },
    ] },
    Element { name: "Rubidium", symbol: "Rb", atomic_number: 37, isotopes: &[
        Isotope { mass_number: 85, mass: 84.911789737999996, abundance: 0.7217 },
    ] },
    Element { name: "Strontium", symbol: "Sr", atomic_number: 38, isotopes: &[
        Isotope { mass_number: 84, mass: 83.913425000000004, abundance: 0.005600000000000001 },
        Isotope { mass_number: 86, mass: 85.909260730900002, abundance: 0.0986 },
        Isotope { mass_number: 87, mass: 86.908877497000006, abundance: 0.07 },
        Isotope { mass_number: 88, mass: 87.905612257100003, abundance: 0.8258 },
    ] },
    Element { name: "Yttrium", symbol: "Y", atomic_number: 39, isotopes: &[
        Isotope { mass_number: 89, mass: 88.905850000000001, abundance: 1.0 },
    ] },
    Element { name: "Zirconium", symbol: "Zr", atomic_number: 40, isotopes: &[
        Isotope { mass_number: 90, mass: 89.9047044, abundance: 0.5145000000000001 },
        Isotope { mass_number: 91, mass: 90.905645800000002, abundance: 0.11220000000000001 },
        Isotope { mass_number: 92, mass: 91.905040799999995, abundance: 0.17149999999999999 },
        Isotope { mass_number: 94, mass: 93.906315199999995, abundance: 0.17379999999999998 },
        Isotope { mass_number: 96, mass: 95.9082776, abundance: 0.0280 },
    ] },
    Element { name: "Nibium", symbol: "Nb", atomic_number: 41, isotopes: &[
        Isotope { mass_number: 93, mass: 92.906378099999998, abundance: 1.0 },
    ] },
    Element { name: "Molybdenum", symbol: "Mo", atomic_number: 42, isotopes: &[
        Isotope { mass_number: 92, mass: 91.906809999999993, abundance: 0.1484 },
        Isotope { mass_number: 94, mass: 93.905088000000006, abundance: 0.0925 },
        Isotope { mass_number: 95, mass: 94.905840999999995, abundance: 0.1592 },
        Isotope { mass_number: 96, mass: 95.904679000000002, abundance: 0.1668 },
        Isotope { mass_number: 97, mass: 96.906020999999996, abundance: 0.0955 },
        Isotope { mass_number: 98, mass: 97.905407999999994, abundance: 0.2413 },
        Isotope { mass_number: 100, mass: 99.907477, abundance: 0.09630000000000001 },
    ] },
    Element { name: "Ruthenium", symbol: "Ru", atomic_number: 44, isotopes: &[
        Isotope { mass_number: 96, mass: 95.907597999999993, abundance: 0.0554 },
        Isotope { mass_number: 98, mass: 97.905287000000001, abundance: 0.0187 },
        Isotope { mass_number: 99, mass: 98.9059393, abundance: 0.1276 },
        Isotope { mass_number: 100, mass: 99.904219499999996, abundance: 0.126 },
        Isotope { mass_number: 101, mass: 100.905582100000004, abundance: 0.17059999999999997 },
        Isotope { mass_number: 102, mass: 101.904349300000007, abundance: 0.3155 },
        Isotope { mass_number: 104, mass: 103.905433000000002, abundance: 0.1862 },
    ] },
    Element { name: "Rhodium", symbol: "Rh", atomic_number: 45, isotopes: &[
        Isotope { mass_number: 103, mass: 102.905500000000004, abundance: 1.0 },
    ] },
    Element { name: "Palladium", symbol: "Pd", atomic_number: 46, isotopes: &[
        Isotope { mass_number: 102, mass: 101.905608999999998, abundance: 0.0102 },
        Isotope { mass_number: 104, mass: 103.904036000000005, abundance: 0.1114 },
        Isotope { mass_number: 105, mass: 104.905085, abundance: 0.22329999999999997 },
        Isotope { mass_number: 106, mass: 105.903486000000001, abundance: 0.2733 },
        Isotope { mass_number: 108, mass: 107.903891999999999, abundance: 0.2646 },
        Isotope { mass_number: 110, mass: 109.905152999999999, abundance: 0.11720000000000001 },
    ] },
    Element { name: "Silver", symbol: "Ag", atomic_number: 47, isotopes: &[
        Isotope { mass_number: 107, mass: 106.905092999999994, abundance: 0.51839 },
        Isotope { mass_number: 109, mass: 108.904756000000006, abundance: 0.48161000000000004 },
    ] },
    Element { name: "Cadmium", symbol: "Cd", atomic_number: 48, isotopes: &[
        Isotope { mass_number: 106, mass: 105.906458000000001, abundance: 0.0125 },
        Isotope { mass_number: 108, mass: 107.904184000000001, abundance: 0.0089 },
        Isotope { mass_number: 110, mass: 109.903002099999995, abundance: 0.1249 },
        Isotope { mass_number: 111, mass: 110.904178099999996, abundance: 0.128 },
        Isotope { mass_number: 112, mass: 111.902757800000003, abundance: 0.2413 },
        Isotope { mass_number: 113, mass: 112.904401699999994, abundance: 0.1222 },
        Isotope { mass_number: 114, mass: 113.903358499999996, abundance: 0.2873 },
        Isotope { mass_number: 116, mass: 115.904756000000006, abundance: 0.07490000000000001 },
    ] },
    Element { name: "Indium", symbol: "In", atomic_number: 49, isotopes: &[
        Isotope { mass_number: 113, mass: 112.904060000000001, abundance: 0.0429 },
        Isotope { mass_number: 115, mass: 114.903878000000006, abundance: 0.9571 },
    ] },
    Element { name: "Tin", symbol: "Sn", atomic_number: 50, isotopes: &[
        Isotope { mass_number: 112, mass: 111.904818000000006, abundance: 0.0097 },
        Isotope { mass_number: 114, mass: 113.902777900000004, abundance: 0.0066 },
        Isotope { mass_number: 115, mass: 114.903341999999995, abundance: 0.0034000000000000002 },
        Isotope { mass_number: 116, mass: 115.901741000000001, abundance: 0.1454 },
        Isotope { mass_number: 117, mass: 116.902951999999999, abundance: 0.0768 },
        Isotope { mass_number: 118, mass: 117.901602999999994, abundance: 0.2422 },
        Isotope { mass_number: 119, mass: 118.903307999999996, abundance: 0.0859 },
        Isotope { mass_number: 120, mass: 119.902194699999996, abundance: 0.3258 },
        Isotope { mass_number: 122, mass: 121.903439000000006, abundance: 0.0463 },
        Isotope { mass_number: 124, mass: 123.905273899999997, abundance: 0.0579 },
    ] },
    Element { name: "Antimony", symbol: "Sb", atomic_number: 51, isotopes: &[
        Isotope { mass_number: 121, mass: 120.903815699999996, abundance: 0.5721 },
        Isotope { mass_number: 123, mass: 122.904213999999996, abundance: 0.4279 },
    ] },
    Element { name: "Tellurium", symbol: "Te", atomic_number: 52, isotopes: &[
        Isotope { mass_number: 120, mass: 119.904020000000003, abundance: 0.0009 },
        Isotope { mass_number: 122, mass: 121.9030439, abundance: 0.0255 },
        Isotope { mass_number: 124, mass: 123.902817900000002, abundance: 0.047400000000000005 },
        Isotope { mass_number: 125, mass: 124.904430700000006, abundance: 0.0707 },
        Isotope { mass_number: 126, mass: 125.903311700000003, abundance: 0.1884 },
        Isotope { mass_number: 128, mass: 127.904463100000001, abundance: 0.31739999999999996 },
        Isotope { mass_number: 130, mass: 129.906224400000014, abundance: 0.3408 },
    ] },
    Element { name: "Iodine", symbol: "I", atomic_number: 53, isotopes: &[
        Isotope { mass_number: 127, mass: 126.904472999999996, abundance: 1.0 },
    ] },
    Element { name: "Xenon", symbol: "Xe", atomic_number: 54, isotopes: &[
        Isotope { mass_number: 128, mass: 127.903531000000001, abundance: 0.0191 },
        Isotope { mass_number: 129, mass: 128.904779999999988, abundance: 0.264 },
        Isotope { mass_number: 130, mass: 129.903509000000014, abundance: 0.040999999999999995 },
        Isotope { mass_number: 131, mass: 130.90507199999999, abundance: 0.212 },
        Isotope { mass_number: 132, mass: 131.904144000000002, abundance: 0.26899999999999996 },
        Isotope { mass_number: 134, mass: 133.905394999999999, abundance: 0.10400000000000001 },
        Isotope { mass_number: 136, mass: 135.90721400000001, abundance: 0.08900000000000001 },
    ] },
    Element { name: "Caesium", symbol: "Cs", atomic_number: 55, isotopes: &[
        Isotope { mass_number: 133, mass: 132.905451932999995, abundance: 1.0 },
    ] },
    Element { name: "Barium", symbol: "Ba", atomic_number: 56, isotopes: &[
        Isotope { mass_number: 132, mass: 131.9050613, abundance: 0.00101 },
        Isotope { mass_number: 134, mass: 133.904508399999997, abundance: 0.024169999999999997 },
        Isotope { mass_number: 135, mass: 134.905688599999991, abundance: 0.06591999999999999 },
        Isotope { mass_number: 136, mass: 135.904575899999998, abundance: 0.07854 },
        Isotope { mass_number: 137, mass: 136.905827399999993, abundance: 0.11231999999999999 },
        Isotope { mass_number: 138, mass: 137.905247199999991, abundance: 0.71698 },
    ] },
    Element { name: "Lanthanum", symbol: "La", atomic_number: 57, isotopes: &[
        Isotope { mass_number: 138, mass: 137.907112000000012, abundance: 0.00089 },
        Isotope { mass_number: 139, mass: 138.906353300000006, abundance: 0.99911 },
    ] },
    Element { name: "Cerium", symbol: "Ce", atomic_number: 58, isotopes: &[
        Isotope { mass_number: 136, mass: 135.907172000000003, abundance: 0.00185 },
        Isotope { mass_number: 138, mass: 137.905991, abundance: 0.00251 },
        Isotope { mass_number: 140, mass: 139.905438699999991, abundance: 0.8845000000000001 },
        Isotope { mass_number: 142, mass: 141.909244000000001, abundance: 0.11114 },
    ] },
    Element { name: "Praseodymium", symbol: "Pr", atomic_number: 59, isotopes: &[
        Isotope { mass_number: 141, mass: 140.907646999999997, abundance: 1.0 },
    ] },
    Element { name: "Neodymium", symbol: "Nd", atomic_number: 60, isotopes: &[
        Isotope { mass_number: 142, mass: 141.907723299999987, abundance: 0.272 },
        Isotope { mass_number: 143, mass: 142.909814299999994, abundance: 0.122 },
        Isotope { mass_number: 144, mass: 143.910087299999987, abundance: 0.23800000000000002 },
        Isotope { mass_number: 145, mass: 144.912573600000002, abundance: 0.083 },
        Isotope { mass_number: 146, mass: 145.913116900000006, abundance: 0.172 },
        Isotope { mass_number: 148, mass: 147.916892999999988, abundance: 0.057999999999999996 },
        Isotope { mass_number: 150, mass: 149.920891000000012, abundance: 0.055999999999999994 },
    ] },
    Element { name: "Samarium", symbol: "Sm", atomic_number: 62, isotopes: &[
        Isotope { mass_number: 144, mass: 143.911999000000009, abundance: 0.0308 },
        Isotope { mass_number: 147, mass: 146.9148979, abundance: 0.15 },
        Isotope { mass_number: 148, mass: 147.914822700000002, abundance: 0.1125 },
        Isotope { mass_number: 149, mass: 148.917184700000007, abundance: 0.1382 },
        Isotope { mass_number: 150, mass: 149.917275499999988, abundance: 0.0737 },
        Isotope { mass_number: 152, mass: 151.919732399999987, abundance: 0.26739999999999997 },
        Isotope { mass_number: 154, mass: 153.92220929999999, abundance: 0.2274 },
    ] },
    Element { name: "Europium", symbol: "Eu", atomic_number: 63, isotopes: &[
        Isotope { mass_number: 151, mass: 150.919857, abundance: 0.4781 },
        Isotope { mass_number: 153, mass: 152.921237, abundance: 0.5219 },
    ] },
    Element { name: "Gadolinium", symbol: "Gd", atomic_number: 64, isotopes: &[
        Isotope { mass_number: 152, mass: 151.919791000000004, abundance: 0.002 },
        Isotope { mass_number: 154, mass: 153.920865600000013, abundance: 0.0218 },
        Isotope { mass_number: 155, mass: 154.92262199999999, abundance: 0.14800000000000002 },
        Isotope { mass_number: 156, mass: 155.922122699999989, abundance: 0.2047 },
        Isotope { mass_number: 157, mass: 156.923960099999988, abundance: 0.1565 },
        Isotope { mass_number: 158, mass: 157.924103900000006, abundance: 0.2484 },
        Isotope { mass_number: 160, mass: 159.927054099999992, abundance: 0.2186 },
    ] },
    Element { name: "Terbium", symbol: "Tb", atomic_number: 65, isotopes: &[
        Isotope { mass_number: 159, mass: 158.925354, abundance: 1.0 },
    ] },
    Element { name: "Dysprosium", symbol: "Dy", atomic_number: 66, isotopes: &[
        Isotope { mass_number: 156, mass: 155.924284, abundance: 0.00056 },
        Isotope { mass_number: 158, mass: 157.92441, abundance: 0.00095 },
        Isotope { mass_number: 160, mass: 159.925203, abundance: 0.02329 },
        Isotope { mass_number: 161, mass: 160.926939, abundance: 0.18889 },
        Isotope { mass_number: 162, mass: 161.926804, abundance: 0.25475 },
        Isotope { mass_number: 163, mass: 162.928737, abundance: 0.24896 },
        Isotope { mass_number: 164, mass: 163.929181, abundance: 0.28260 },
    ] },
    Element { name: "Holmium", symbol: "Ho", atomic_number: 67, isotopes: &[
        Isotope { mass_number: 165, mass: 164.930328, abundance: 1.0 },
    ] },
    Element { name: "Erbium", symbol: "Er", atomic_number: 68, isotopes: &[
        Isotope { mass_number: 162, mass: 161.928787, abundance: 0.00056 },
        Isotope { mass_number: 164, mass: 163.929207, abundance: 0.01601 },
        Isotope { mass_number: 166, mass: 165.930299, abundance: 0.33503 },
        Isotope { mass_number: 167, mass: 166.932054, abundance: 0.22869 },
        Isotope { mass_number: 168, mass: 167.932376, abundance: 0.26978 },
        Isotope { mass_number: 170, mass: 169.93547, abundance: 0.14910 },
    ] },
    Element { name: "Thulium", symbol: "Tm", atomic_number: 69, isotopes: &[
        Isotope { mass_number: 169, mass: 168.934218, abundance: 1.0 },
    ] },
    Element { name: "Ytterbium", symbol: "Yb", atomic_number: 70, isotopes: &[
        Isotope { mass_number: 168, mass: 167.933889, abundance: 0.00126 },
        Isotope { mass_number: 170, mass: 169.93476725, abundance: 0.03023 },
        Isotope { mass_number: 171, mass: 170.93633152, abundance: 0.14216 },
        Isotope { mass_number: 172, mass: 171.93638666, abundance: 0.21754 },
        Isotope { mass_number: 173, mass: 172.93821622, abundance: 0.16098 },
        Isotope { mass_number: 174, mass: 173.93886755, abundance: 0.31896 },
        Isotope { mass_number: 176, mass: 175.9425747, abundance: 0.12887 },
    ] },
    Element { name: "Lutetium", symbol: "Lu", atomic_number: 71, isotopes: &[
        Isotope { mass_number: 175, mass: 174.940777, abundance: 0.97401 },
        Isotope { mass_number: 176, mass: 175.942692, abundance: 0.02599 },
    ] },
    Element { name: "Hafnium", symbol: "Hf", atomic_number: 72, isotopes: &[
        Isotope { mass_number: 176, mass: 175.941408599999988, abundance: 0.0526 },
        Isotope { mass_number: 177, mass: 176.943220700000012, abundance: 0.18600000000000003 },
        Isotope { mass_number: 178, mass: 177.943698799999993, abundance: 0.2728 },
        Isotope { mass_number: 179, mass: 178.945816100000002, abundance: 0.1362 },
        Isotope { mass_number: 180, mass: 179.946550000000002, abundance: 0.3508 },
    ] },
    Element { name: "Tantalum", symbol: "Ta", atomic_number: 73, isotopes: &[
        Isotope { mass_number: 180, mass: 179.94747, abundance: 0.0001176 },
        Isotope { mass_number: 181, mass: 180.947995800000001, abundance: 0.99988 },
    ] },
    Element { name: "Tungsten", symbol: "W", atomic_number: 74, isotopes: &[
        Isotope { mass_number: 180, mass: 179.946704000000011, abundance: 0.0012 },
        Isotope { mass_number: 182, mass: 181.948204199999992, abundance: 0.265 },
        Isotope { mass_number: 183, mass: 182.950222999999994, abundance: 0.1431 },
        Isotope { mass_number: 184, mass: 183.950930999999997, abundance: 0.3064 },
        Isotope { mass_number: 186, mass: 185.954364099999992, abundance: 0.2843 },
    ] },
    Element { name: "Rhenium", symbol: "Re", atomic_number: 75, isotopes: &[
        Isotope { mass_number: 185, mass: 184.952955000000003, abundance: 0.374 },
        Isotope { mass_number: 187, mass: 186.95575310000001, abundance: 0.626 },
    ] },
    Element { name: "Osmium", symbol: "Os", atomic_number: 76, isotopes: &[
        Isotope { mass_number: 184, mass: 183.952493, abundance: 0.0002 },
        Isotope { mass_number: 186, mass: 185.953838, abundance: 0.0159 },
        Isotope { mass_number: 187, mass: 186.955750, abundance: 0.0196 },
        Isotope { mass_number: 188, mass: 187.955837, abundance: 0.1324 },
        Isotope { mass_number: 189, mass: 188.958146, abundance: 0.1615 },
        Isotope { mass_number: 190, mass: 189.958446, abundance: 0.2626 },
        Isotope { mass_number: 192, mass: 191.96148, abundance: 0.4078 },
    ] },
    Element { name: "Iridium", symbol: "Ir", atomic_number: 77, isotopes: &[
        Isotope { mass_number: 191, mass: 190.960591, abundance: 0.3723 },
        Isotope { mass_number: 193, mass: 192.962924, abundance: 0.6277 },
    ] },
    Element { name: "Platinum", symbol: "Pt", atomic_number: 78, isotopes: &[
        Isotope { mass_number: 192, mass: 191.961038000000002, abundance: 0.00782 },
        Isotope { mass_number: 194, mass: 193.962680299999988, abundance: 0.32966999999999996 },
        Isotope { mass_number: 195, mass: 194.964791100000014, abundance: 0.33832 },
        Isotope { mass_number: 196, mass: 195.964951500000012, abundance: 0.25242000000000003 },
        Isotope { mass_number: 198, mass: 197.967893000000004, abundance: 0.07163 },
    ] },
    Element { name: "Gold", symbol: "Au", atomic_number: 79, isotopes: &[
        Isotope { mass_number: 197, mass: 196.96655100000001, abundance: 1.0 },
    ] },
    Element { name: "Mercury", symbol: "Hg", atomic_number: 80, isotopes: &[
        Isotope { mass_number: 196, mass: 195.965833000000004, abundance: 0.0015 },
        Isotope { mass_number: 198, mass: 197.966768999999999, abundance: 0.09970000000000001 },
        Isotope { mass_number: 199, mass: 198.968279899999999, abundance: 0.16870000000000002 },
        Isotope { mass_number: 200, mass: 199.968325999999991, abundance: 0.231 },
        Isotope { mass_number: 201, mass: 200.970302299999986, abundance: 0.1318 },
        Isotope { mass_number: 202, mass: 201.970642999999996, abundance: 0.2986 },
        Isotope { mass_number: 204, mass: 203.973493899999994, abundance: 0.0687 },
    ] },
    Element { name: "Thallium", symbol: "Tl", atomic_number: 81, isotopes: &[
        Isotope { mass_number: 203, mass: 202.972344200000009, abundance: 0.2952 },
        Isotope { mass_number: 205, mass: 204.97442749999999, abundance: 0.7048000000000001 },
    ] },
    Element { name: "Lead", symbol: "Pb", atomic_number: 82, isotopes: &[
        Isotope { mass_number: 204, mass: 203.973043600000011, abundance: 0.013999999999999999 },
        Isotope { mass_number: 206, mass: 205.974465299999991, abundance: 0.24100000000000002 },
        Isotope { mass_number: 207, mass: 206.975896900000009, abundance: 0.221 },
        Isotope { mass_number: 208, mass: 207.976653800000008, abundance: 0.524 },
    ] },
    Element { name: "Bismuth", symbol: "Bi", atomic_number: 83, isotopes: &[
        Isotope { mass_number: 209, mass: 208.980398699999995, abundance: 1.0 },
    ] },
    Element { name: "Thorium", symbol: "Th", atomic_number: 90, isotopes: &[
        Isotope { mass_number: 230, mass: 230.033133800000002, abundance: 0.0002 },
        Isotope { mass_number: 232, mass: 232.038055299999996, abundance: 0.9998 },
    ] },
    Element { name: "Protactinium", symbol: "Pa", atomic_number: 91, isotopes: &[
        Isotope { mass_number: 231, mass: 231.03588, abundance: 1.0 },
    ] },
    Element { name: "Uranium", symbol: "U", atomic_number: 92, isotopes: &[
        Isotope { mass_number: 234, mass: 234.040950, abundance: 0.000054 },
        Isotope { mass_number: 235, mass: 235.043928, abundance: 0.007204 },
        Isotope { mass_number: 238, mass: 238.05079, abundance: 0.992742 },
    ] },
];
