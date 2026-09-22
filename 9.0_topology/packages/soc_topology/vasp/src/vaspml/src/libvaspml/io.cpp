#include "io.hpp"
#include "SemanticVersion.hpp"

using namespace vaspml;

void io::checkLine(const String& line, const String& search, const Int n)
{
    if (!string_tools::checkForSubstring(line, search))
    {
        throw std::runtime_error("ERROR: Invalid file, line " + std::to_string(n)
                                 + ": expected line containing string\n   \"" + search
                                 + "\"\nbut got instead:\n   \"" + line + "\"\n");
    }

    return;
}

Int io::readLine(std::ifstream& input, std::string& buffer, Int skip, bool expectEOF)
{
    for (Int i = 0; i < skip + 1; ++i)
    {
        if (!getline(input, buffer))
        {
            if (expectEOF) return i;
            else throw std::runtime_error("ERROR: Unexpected end of file.");
        }
    }

    return skip + 1;
}

void io::writeJson(const Record& record, std::ostream& output, Int level, Int indentSpaces)
{
    const String     indent1 = String(indentSpaces, ' ');
    const String     indent = String(indentSpaces * (level + 1), ' ');
    const Vec1String keys = record.keys();
    output << "{\n";
    const auto last = keys.cend() - 1;
    for (auto k = keys.cbegin(); k != keys.cend(); ++k)
    {
        output << indent + '"' + *k + "\" : ";
        std::visit(
            [&](auto&& arg)
            {
                using T = std::decay_t<decltype(arg)>;
                //************************************
                // Item types: Real, Int, String, Bool
                //************************************
                if constexpr (isScalar<T>::value) writeJsonScalarItem<T>(output, arg);
                //******************
                // Item types: ShRec
                //******************
                else if constexpr (std::is_same_v<T, ShRec>)
                {
                    io::writeJson(*arg, output, level + 1, indentSpaces);
                    output << indent << '}';
                }
                else if constexpr (isVector<T>::value)
                {
                    output << "[\n";
                    using T1 = typename T::value_type;
                    auto lasta1 = arg.cend() - 1;
                    for (auto a1 = arg.cbegin(); a1 != arg.cend(); ++a1)
                    {
                        output << indent << indent1;
                        //******************************************
                        // Item types: Vec1Real, Vec1Int, Vec1String
                        //******************************************
                        if constexpr (isScalar<T1>::value) writeJsonScalarItem<T1>(output, *a1);
                        //**********************
                        // Item types: Vec1ShRec
                        //**********************
                        else if constexpr (std::is_same_v<T1, ShRec>)
                        {
                            io::writeJson(**a1, output, level + 2, indentSpaces);
                            output << indent << indent1 << '}';
                        }
                        //******************************
                        // Item types: Vec2Real, Vec2Int
                        //******************************
                        else if constexpr (isVector<T1>::value)
                        {
                            output << "[\n";
                            using T2 = typename T1::value_type;
                            auto lasta2 = a1->cend() - 1;
                            for (auto a2 = a1->cbegin(); a2 != a1->cend(); ++a2)
                            {
                                output << indent << indent1 << indent1;
                                writeJsonScalarItem<T2>(output, *a2);
                                if (a2 != lasta2) output << ",\n";
                            }
                            output << '\n' << indent << indent1 << ']';
                        }
                        else output << "null";
                        if (a1 != lasta1) output << ",\n";
                    }
                    output << '\n' << indent << ']';
                }
                else output << "null";
                if (k != last) output << ',';
                output << '\n';
            },
            record.at(*k));
    }
    if (level == 0) output << "}\n";

    return;
}

void io::readMlab(Record& record, std::ifstream& input)
{
    using namespace string_tools;
    String line;
    Int    n = 0;

    // Read first line and check for version string.
    n += readLine(input, line);
    if (!checkForSubstring(line, "Version"))
    {
        std::runtime_error("ERROR: Invalid header found in input file in ML_AB format.");
    }
    // Store version in record.
    Vec1String      split = splitString(trim(line), " ");
    SemanticVersion version(split.at(0));
    record["version"] = version.toString();

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The number of configurations", n););
    n += readLine(input, line, 1);
    record["numStructures"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atom type", n););
    n += readLine(input, line, 1);
    record["maxTypes"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The atom types in the data file", n););
    n += readLine(input, line, 0);
    record["types"] = readMlabVector<Vec1String>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atoms per system", n););
    n += readLine(input, line, 1);
    record["maxAtoms"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atoms per atom type", n););
    n += readLine(input, line, 1);
    record["maxAtomsPerType"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "Reference atomic energy (eV)", n););
    n += readLine(input, line, 0);
    record["atomicRefEnergy"] = readMlabVector<Vec1Real>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "Atomic mass", n););
    n += readLine(input, line, 0);
    record["atomicMass"] = readMlabVector<Vec1Real>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "The numbers of basis sets per atom type", n););
    n += readLine(input, line, 0);
    record["numLrc"] = readMlabVector<Vec1Int>(input, n);
    const Int maxTypes = record.cget<Int>("maxTypes");
    record["lrcStructure"] = Vec2Int(maxTypes);
    record["lrcAtom"] = Vec2Int(maxTypes);

    for (Int i = 0; i < maxTypes; ++i)
    {
        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Basis set for", n););
        n += readLine(input, line, 0);
        const Int numLrc = record.cget<Vec1Int>("numLrc").at(i);
        Vec1Int   lrcStructure;
        Vec1Int   lrcAtom;
        for (Int j = 0; j < numLrc; ++j)
        {
            n += readLine(input, line, 0);
            Vec1String values = splitString(line, " ");
            lrcStructure.push_back(readMlabScalar<Int>(values.at(0)));
            lrcAtom.push_back(readMlabScalar<Int>(values.at(1)));
        }
        record.get<Vec2Int>("lrcStructure").at(i) = lrcStructure;
        record.get<Vec2Int>("lrcAtom").at(i) = lrcAtom;
        n += readLine(input, line, 0);
    }

    record["structures"] = Vec1ShRec();
    Vec1ShRec& structures = record.get<Vec1ShRec>("structures");
    const Int  numStructures = record.cget<Int>("numStructures");
    for (Int i = 0; i < numStructures; ++i)
    {
        structures.push_back(std::make_shared<Record>());
        Record& structure = *structures.back();
        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Configuration num.", n););

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "System name", n););
        n += readLine(input, line, 1);
        structure["system"] = readMlabScalar<String>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "The number of atom types", n););
        n += readLine(input, line, 1);
        structure["numTypes"] = readMlabScalar<Int>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "The number of atoms", n););
        n += readLine(input, line, 1);
        structure["numAtoms"] = readMlabScalar<Int>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "Atom types and atom numbers", n););
        n += readLine(input, line, 0);
        const Int  numTypes = structure.cget<Int>("numTypes");
        Vec1String types;
        Vec1Int    numAtomsPerType;
        for (Int j = 0; j < numTypes; ++j)
        {
            n += readLine(input, line, 0);
            Vec1String values = splitString(line, " ");
            types.push_back(readMlabScalar<String>(values.at(0)));
            numAtomsPerType.push_back(readMlabScalar<Int>(values.at(1)));
        }
        structure["types"] = types;
        structure["numAtomsPerType"] = numAtomsPerType;

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "CTIFOR", n););
        n += readLine(input, line, 1);
        structure["CTIFOR"] = readMlabScalar<Real>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "Primitive lattice vectors (ang.)", n););
        n += readLine(input, line, 0);
        structure["lattice"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Atomic positions (ang.)", n););
        n += readLine(input, line, 0);
        structure["positions"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Total energy (eV)", n););
        n += readLine(input, line, 1);
        structure["energy"] = readMlabScalar<Real>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "Forces (eV ang.^-1)", n););
        n += readLine(input, line, 0);
        structure["forces"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Stress (kbar)", n););
        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "XX YY ZZ", n););
        n += readLine(input, line, 0);
        Vec1Real stress = readMlabVector<Vec1Real>(input, n);
        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "XY YZ ZX", n););
        n += readLine(input, line, 0);
        Vec1Real stress2 = readMlabVector<Vec1Real>(input, n, (i == numStructures - 1));
        stress.insert(stress.end(), stress2.begin(), stress2.end());
        structure["stress"] = stress;
    }

    return;
}
