/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_AUTOGENERATEDHELPER_H
#define LEVI_AUTOGENERATEDHELPER_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/TypeDetector.h>
#include <levi/TreeExpander.h>
#include <levi/Expression.h>
#include <levi/autogenerated/Path.h>

#include <levi/external/zupply.h>
#include <shlibpp/SharedLibraryClass.h>
#include <shlibpp/SharedLibrary.h>

#include <ostream>
#include <cstdlib>

#include <fstream>

//Taken from https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
class static_string
{
    const char* const p_;
    const std::size_t sz_;

public:

    template <std::size_t N>
    constexpr static_string(const char(&a)[N]) noexcept
        : p_(a)
          , sz_(N-1)
    {}

    constexpr static_string(const char* p, std::size_t N) noexcept
        : p_(p)
          , sz_(N)
    {}

    constexpr const char* data() const noexcept {return p_;}
    constexpr std::size_t size() const noexcept {return sz_;}

};

inline std::ostream& operator<<(std::ostream& os, static_string const& s)
{
    return os.write(s.data(), static_cast<std::streamsize>(s.size()));
}

template <class T>
constexpr static_string type_name()
{
#ifdef __clang__
    static_string p = __PRETTY_FUNCTION__;
    return static_string(p.data() + 31, p.size() - 31 - 1);
#elif defined(__GNUC__)
    static_string p = __PRETTY_FUNCTION__;
#  if __cplusplus < 201402
    return static_string(p.data() + 36, p.size() - 36 - 1);
#  else
    return static_string(p.data() + 46, p.size() - 46 - 1);
#  endif
#elif defined(_MSC_VER)
    static_string p = __FUNCSIG__;
    return static_string(p.data() + 38, p.size() - 38 - 7);
#endif
}

namespace levi {
    template <typename EvaluableT>
    class AutogeneratedHelper;

    template <typename BaseType>
    class CompiledEvaluableFactory;
}

template <typename BaseType>
class levi::CompiledEvaluableFactory {

    template <typename EvaluableT>
    friend class AutogeneratedHelper;

    BaseType* m_compiledEvaluable;
    shlibpp::SharedLibraryClassFactory<BaseType> m_compiledEvaluableFactory;

public:

    CompiledEvaluableFactory()
        : m_compiledEvaluable(nullptr)
    {}

    ~CompiledEvaluableFactory() {
        if (m_compiledEvaluable) {
            m_compiledEvaluableFactory.destroy(m_compiledEvaluable);
            m_compiledEvaluable = nullptr;
        }
    }

    BaseType* operator->() const {
        return m_compiledEvaluable;
    }
};

template <typename EvaluableT>
class levi::AutogeneratedHelper {
public:

    using SqueezedMatrix = typename levi::TreeComponent<EvaluableT>::SqueezedMatrix;
    using SqueezedMatrixRef = Eigen::Ref<SqueezedMatrix>;

private:

    using Type = levi::EvaluableType;

    struct LiteralComponent {
        std::string literal;
        bool isScalar = false;
        bool isAlreadyCompressed = false;
        bool isSimple = false;
    };

    struct CommonSubExp {
        std::string name;
        std::string expression;
        std::string declaration;
        bool active = true;
    };

    std::string m_cleanName;
    std::string m_workingDirectory;
    std::vector<LiteralComponent> m_literalSubExpressions;
    std::vector<levi::TreeComponent<EvaluableT>> m_expandedExpression;
    std::vector<size_t> m_generics;
    std::vector<size_t> m_finalExpressionIndices;
    std::unordered_map<std::string, std::vector<int>> m_helpersMap;
    std::vector<CommonSubExp> m_helpers;
    std::ostringstream m_helpersDeclarations;
    std::vector<CommonSubExp> m_commonExpressions;
    std::ostringstream m_commonDeclarations;
    std::vector<std::ostringstream> m_finalExpressions;
    std::string m_genericsName, m_helpersName, m_commonsName;
    std::unordered_map<std::string, std::vector<int>> m_commonsMap;

    std::vector<SqueezedMatrixRef> m_genericsRefs;

    void expandElement(int i) {
        Type type;

        levi::TreeComponent<EvaluableT>& subExpr = m_expandedExpression[static_cast<size_t>(i)];
        LiteralComponent& literalSubExpr = m_literalSubExpressions[static_cast<size_t>(i)];
        const LiteralComponent& lhs = m_literalSubExpressions[subExpr.lhsIndex];
        const LiteralComponent& rhs = m_literalSubExpressions[subExpr.rhsIndex];

        if (literalSubExpr.isAlreadyCompressed) {
            return;
        }

        type = subExpr.type;

        if (type == Type::Sum) {

            literalSubExpr.literal = "(" + lhs.literal + " + " + rhs.literal + ")";

            literalSubExpr.isScalar = (lhs.isScalar && rhs.isScalar);
            literalSubExpr.isSimple = false;

        } else if (type == Type::Subtraction) {

            literalSubExpr.literal = "(" + lhs.literal + " - " + rhs.literal + ")";
            literalSubExpr.isScalar = (lhs.isScalar && rhs.isScalar);
            literalSubExpr.isSimple = false;


        } else if (type == Type::Product) {

            literalSubExpr.literal = lhs.literal + " * " + rhs.literal;
            literalSubExpr.isScalar = (lhs.isScalar && rhs.isScalar);
            literalSubExpr.isSimple = false;


        } else if (type == Type::Division) {

            literalSubExpr.literal = lhs.literal + " / " + rhs.literal;
            literalSubExpr.isScalar = lhs.isScalar;
            literalSubExpr.isSimple = false;


        } else if (type == Type::InvertedSign) {

            literalSubExpr.literal = "-" + lhs.literal;
            literalSubExpr.isScalar = lhs.isScalar;
            literalSubExpr.isSimple = false;


        } else if (type == Type::Pow) {

            literalSubExpr.literal = "std::pow(" + lhs.literal +", " + std::to_string(subExpr.exponent) + ")";
            literalSubExpr.isScalar = true;
            literalSubExpr.isSimple = false;


        } else if (type == Type::Transpose) {

            assert(!lhs.isScalar);
            literalSubExpr.literal = "(" + lhs.literal + ").transpose()";
            literalSubExpr.isScalar = false;
            literalSubExpr.isSimple = true;


        } else if (type == Type::Row) {

            assert(!lhs.isScalar);
            if ((subExpr.rows() == 1) && (subExpr.cols() == 1)) {
                subExpr.block.startCol = 0;
                literalSubExpr.literal = "(" + lhs.literal + ")(" + std::to_string(subExpr.block.startRow) + ", 0)";
                subExpr.type = Type::Element;
                literalSubExpr.isScalar = true;
                literalSubExpr.isSimple = true;
            } else {
                literalSubExpr.literal = "(" + lhs.literal + ").row(" + std::to_string(subExpr.block.startRow) + ")";
                literalSubExpr.isScalar = false;
                literalSubExpr.isSimple = true;
            }
        } else if (type == Type::Column) {

            assert(!lhs.isScalar);
            if ((subExpr.rows() == 1) && (subExpr.cols() == 1)) {
                subExpr.block.startRow = 0;
                literalSubExpr.literal = "(" + lhs.literal + ")(0, " + std::to_string(subExpr.block.startCol) + ")";
                subExpr.type = Type::Element;
                literalSubExpr.isScalar = true;
                literalSubExpr.isSimple = true;
            } else {
                literalSubExpr.literal = "(" + lhs.literal + ").col(" + std::to_string(subExpr.block.startCol) + ")";
                literalSubExpr.isScalar = false;
                literalSubExpr.isSimple = true;
            }

        } else if (type == Type::Element) {

            assert(!lhs.isScalar);
            literalSubExpr.literal = "(" + lhs.literal + ")(" + std::to_string(subExpr.block.startRow) + ", " + std::to_string(subExpr.block.startCol) + ")";
            literalSubExpr.isScalar = true;
            literalSubExpr.isSimple = true;


        } else if (type == Type::Block) {

            assert(!lhs.isScalar);
            literalSubExpr.literal = "(" + lhs.literal + ").block<" + std::to_string(subExpr.block.rows) + ", " + std::to_string(subExpr.block.cols) + ">(" + std::to_string(subExpr.block.startRow) + ", " + std::to_string(subExpr.block.startCol) + ")";
            literalSubExpr.isScalar = false;
            literalSubExpr.isSimple = true;


        }

        if (!literalSubExpr.isScalar && (subExpr.rows() == 1) && (subExpr.cols() == 1)) { //the subexpression is an Eigen object but of dimension 1x1
            literalSubExpr.literal = "(" + literalSubExpr.literal + ")(0,0)";
            literalSubExpr.isScalar = true;
        }
    }


    void getScalarVariables() {

        for(int index = 0; index < m_expandedExpression.size(); ++index) {

            if (!m_literalSubExpressions[index].isScalar || m_literalSubExpressions[index].isAlreadyCompressed ||
                m_literalSubExpressions[index].isSimple) {
                continue;
            }

            std::vector<int>& helpers = m_helpersMap[m_literalSubExpressions[index].literal];

            if (helpers.size() > 1) {
                CommonSubExp newHelper;
                newHelper.name = m_helpersName + std::to_string(m_helpers.size()) + "_";
                newHelper.expression = m_literalSubExpressions[index].literal;
                std::ostringstream decl;
                decl << "    " << type_name<typename EvaluableT::value_type>() << " " << newHelper.name << " = ";
                newHelper.declaration = decl.str();

                for (int& i : helpers) {
                    m_literalSubExpressions[i].literal = newHelper.name;
                    m_literalSubExpressions[i].isAlreadyCompressed = true;
                }

                m_helpers.push_back(newHelper);
            }
        }
    }

    void getCommonSubExpression() {

        for(int index = 0; index < m_expandedExpression.size(); ++index) {

            if (m_literalSubExpressions[index].isAlreadyCompressed || m_literalSubExpressions[index].isSimple) {
                continue;
            }

            std::vector<int>& commons = m_commonsMap[m_literalSubExpressions[index].literal];

            if (commons.size() > 1) {
                CommonSubExp newCommon;
                newCommon.name = m_commonsName + std::to_string(m_commonExpressions.size()) + "_";
                newCommon.expression = m_literalSubExpressions[index].literal;
                std::ostringstream decl;
                decl << "    Eigen::Matrix<" << type_name<typename EvaluableT::value_type>() << ", "
                     << std::to_string(m_expandedExpression[index].rows()) << ", "
                     << std::to_string(m_expandedExpression[index].cols()) << "> "
                     << newCommon.name << " = ";
                newCommon.declaration = decl.str();

                for (int& i : commons) {
                    m_literalSubExpressions[i].literal = newCommon.name;
                    m_literalSubExpressions[i].isAlreadyCompressed = true;
                }

                m_commonExpressions.push_back(newCommon);
            }
        }
    }

    void cleanHelpers() {

        for (int i = m_helpers.size() - 1; i >= 0; --i) {

            bool isInFinal = false;
            for (const size_t& finals : m_finalExpressionIndices) {
                std::string& final = m_literalSubExpressions[finals].literal;
                isInFinal = isInFinal || final.find(m_helpers[i].name) != std::string::npos;
            }

            std::vector<std::pair<int, size_t>> inclusions;

            for (int other = i + 1; other < m_helpers.size(); ++other) { //Check if the current helper was present in a subexpression previously checked
                if (m_helpers[other].active && (m_helpers[other].expression.size() > m_helpers[i].expression.size())) {
                    size_t offset = m_helpers[other].expression.find(m_helpers[i].expression);
                    if (offset != std::string::npos) {
                        inclusions.emplace_back(other, offset);
                    }
                }
            }

            if (isInFinal) {
                for (auto& other: inclusions) {
                    m_helpers[other.first].expression.replace(other.second, m_helpers[i].expression.size(), m_helpers[i].name);
                }
            } else {
                if (inclusions.size() > 1) {
                    for (auto& other: inclusions) {
                        m_helpers[other.first].expression.replace(other.second, m_helpers[i].expression.size(), m_helpers[i].name);
                    }
                } else {
                    for (int deactivatedHelper : m_helpersMap[m_helpers[i].expression]) {
                        m_literalSubExpressions[deactivatedHelper].literal = m_helpers[i].expression; //restore the original expression in the deleted helper
                    }
                    m_helpers[i].active = false;
                }
            }
        }

        for (auto& helper :m_helpers) {
            if (helper.active) {
                m_helpersDeclarations << helper.declaration << helper.expression << ";" << std::endl;
            }
        }

    }

    void cleanCommonSubExpressions() {

        for (int i = m_commonExpressions.size() - 1; i >= 0; --i) {

            bool isInFinal = false;
            for (const size_t& finals : m_finalExpressionIndices) {
                std::string& final = m_literalSubExpressions[finals].literal;
                isInFinal = isInFinal || final.find(m_commonExpressions[i].name) != std::string::npos;
            }

            std::vector<std::pair<int, size_t>> inclusions;

            for (int other = i + 1; other < m_commonExpressions.size(); ++other) { //Check if the current subexpression was present in a subexpression previously checked
                if (m_commonExpressions[other].active && (m_commonExpressions[other].expression.size() > m_commonExpressions[i].expression.size())) {
                    size_t offset = m_commonExpressions[other].expression.find(m_commonExpressions[i].expression);
                    if (offset != std::string::npos) {
                        inclusions.emplace_back(other, offset);
                    }
                }
            }

            if (isInFinal) {
                for (auto& other: inclusions) {
                    m_commonExpressions[other.first].expression.replace(other.second, m_commonExpressions[i].expression.size(), m_commonExpressions[i].name);
                }
            } else {
                if (inclusions.size() > 1) {
                    for (auto& other: inclusions) {
                        m_commonExpressions[other.first].expression.replace(other.second, m_commonExpressions[i].expression.size(), m_commonExpressions[i].name);
                    }
                } else {
                    m_commonExpressions[i].active = false;
                }
            }
        }

        for (auto& common :m_commonExpressions) {
            if (common.active) {
                m_commonDeclarations << common.declaration << common.expression << ";" << std::endl;
            }
        }

    }


    void splitExpression(const std::string& expression, std::ostringstream& output) {
        size_t offset = 0;
        size_t currentindex = 0;
        int lastStopCharachter = -1;

        while (currentindex < expression.size()) {

            if (expression[currentindex] == '+' || expression[currentindex] == '-' || expression[currentindex] == '*') {
                lastStopCharachter = currentindex;
            }

            if (lastStopCharachter > 0 && (currentindex - offset) > 150) {
                output << expression.substr(offset, (static_cast<size_t>(lastStopCharachter) - offset + 1)) << std::endl << "       ";
                offset = static_cast<size_t>(lastStopCharachter) + 1;
                lastStopCharachter = -1;
            }

            ++currentindex;

        }

        output << expression.substr(offset, currentindex - offset);
    }


    void getLiteralExpression() {
        m_literalSubExpressions.resize(m_expandedExpression.size());
        m_helpersMap.clear();

        std::cout << "Expanding generics.." << std::endl;


        for (size_t generic = 0; generic < m_generics.size(); ++generic) {
            if (m_expandedExpression[m_generics[generic]].rows() == 1 && m_expandedExpression[m_generics[generic]].cols() == 1) {
                m_literalSubExpressions[m_generics[generic]].literal = m_genericsName + "[" + std::to_string(generic) + "](0,0)";
                m_literalSubExpressions[m_generics[generic]].isScalar = true;
                m_literalSubExpressions[m_generics[generic]].isSimple = true;
            } else {
                m_literalSubExpressions[m_generics[generic]].literal = m_genericsName + "[" + std::to_string(generic) + "]";
                m_literalSubExpressions[m_generics[generic]].isScalar = false;
            }
        }

        std::cout << "First tree expansion.." << std::endl;

        for(int i = 0; i < m_expandedExpression.size(); ++i) {

            expandElement(i);

            LiteralComponent& literalSubExpr = m_literalSubExpressions[static_cast<size_t>(i)];

            m_helpersMap[literalSubExpr.literal].push_back(i);
        }

        std::cout << "Getting helpers.." << std::endl;

        getScalarVariables();

        std::cout << "Second tree expansion.." << std::endl;

        for(int i = 0; i < m_expandedExpression.size(); ++i) {
            expandElement(i); //some inner elements may have been modified when extracting the helpers. Here we reconstruct the total expression
        }

        std::cout << "Cleaning helpers.." << std::endl;

        cleanHelpers(); //remove helpers that turned out to be useless;

        std::cout << "Third tree expansion.." << std::endl;

        for(int i = 0; i < m_expandedExpression.size(); ++i) {
            levi::TreeComponent<EvaluableT>& subExpr = m_expandedExpression[static_cast<size_t>(i)];
            LiteralComponent& literalSubExpr = m_literalSubExpressions[static_cast<size_t>(i)];

            expandElement(i); //some inner elements may have been modified when extracting the helpers. Here we reconstruct the total expression

            m_commonsMap[literalSubExpr.literal].push_back(i);
        }

        std::cout << "Getting common subexpressions.." << std::endl;

        getCommonSubExpression(); //check if some subexpressions are duplicated

        std::cout << "Fourth tree expansion.." << std::endl;

        for(int i = 0; i < m_expandedExpression.size(); ++i) {
            expandElement(i); //some inner elements may have been modified when extracting the common subexpressions. Here we reconstruct the total expression
        }

        std::cout << "Cleaning common subexpressions.." << std::endl;

        cleanCommonSubExpressions(); //remove common subexpressions that turned out to be useless;

        for (size_t final = 0; final < m_finalExpressionIndices.size(); ++final) {

            std::cout << "Adding new lines.." << std::endl;

            splitExpression(m_literalSubExpressions[m_finalExpressionIndices[final]].literal, m_finalExpressions[final]); //add some newlines where needed
        }

    }

public:

    AutogeneratedHelper() { }

    AutogeneratedHelper(const std::vector<levi::ExpressionComponent<EvaluableT>>& fullExpressions, const std::string& name)
    {
        m_genericsName = "generics";
        m_helpersName = "m_helper";
        m_commonsName = "m_common";

        setExpressions(fullExpressions, name);
    }

    void setVariablesName(const std::string& genericsName = "generics", const std::string& helpersName = "m_helper",
                          const std::string& commonsName = "m_common") {
        m_genericsName = genericsName;
        m_helpersName = helpersName;
        m_commonsName = commonsName;
    }

    void setExpressions(const std::vector<levi::ExpressionComponent<EvaluableT>>& fullExpressions, const std::string& name) {

        for (const auto& expressions : fullExpressions) {
            std::cout << "Expanding expression.." << std::endl;
            m_finalExpressionIndices.emplace_back(levi::expandTree(expressions, false, m_expandedExpression, m_generics));
        }
        m_finalExpressions.resize(m_finalExpressionIndices.size());
        getLiteralExpression();

        m_cleanName = name;

        for (char& letter : m_cleanName) {
            if (letter == ' ') {
                letter = '_';
            } else if (!(((letter >= 'a') && (letter <= 'z')) || ((letter >= 'A') && (letter <= 'Z')) || ((letter >= '0') && (letter <= '9')))) {
                letter = '-';
            }
        }

        for (size_t generic : m_generics) {
            m_genericsRefs.emplace_back(m_expandedExpression[generic].buffer);
        }
    }

    bool setWorkingDirectory(std::string workingDirectory) {
        zz::fs::Path dir(workingDirectory);
        if (!dir.exist()) {
            bool dirCreated = zz::os::create_directory(workingDirectory);

            if (dirCreated) {
                m_workingDirectory = workingDirectory;
                return true;
            }

            return false;
        }
        m_workingDirectory = workingDirectory;
        return true;
    }

    const std::string& name() {
        return m_cleanName;
    }

    template <typename BaseClass>
    void compile(const std::string& headerContent, const std::string& cppContent,
                 levi::CompiledEvaluableFactory<BaseClass>& baseClassFactory,
                 const std::string& className) {

        if (!m_workingDirectory.size()) {

            std::cout << "Creating directory " << zz::os::current_working_directory() + "/" + m_cleanName << std::endl;

            bool dirCreated = setWorkingDirectory(zz::os::current_working_directory() + "/" + m_cleanName);
            assert(dirCreated);
        }

        std::string leviListDir = LEVI_AUTOGENERATED_DIR;

        std::string CMakeSource = leviListDir + "/CMakeLists.auto";
        std::string CMakeDest = m_workingDirectory + "/CMakeLists.txt";

        std::cout << "Copying CMakeLists.." << std::endl;

        zz::os::copyfile(CMakeSource, CMakeDest);

        std::string headerName = m_workingDirectory + "/source.h";

        std::cout << "Adding header.." << std::endl;

        std::fstream header(headerName.c_str(), std::ios::out | std::ios::trunc);

        assert(header.is_open());
        header << headerContent;
        header.close();

        std::string cppName = m_workingDirectory + "/source.cpp";

        std::cout << "Adding cpp.." << std::endl;

        std::fstream cpp(cppName.c_str(), std::ios::out | std::ios::trunc);
        assert(cpp.is_open());

        cpp << "//This file has been autogenerated" << std::endl;
        cpp << "#include <shlibpp/SharedLibraryClass.h>" << std::endl;
        cpp << "#include \"source.h\" " << std::endl;
        cpp << "typedef " << type_name<BaseClass>() << " base_type;" << std::endl;
        cpp << "SHLIBPP_DEFINE_SHARED_SUBCLASS(" << className << "Factory, "  << className <<", base_type);"  << std::endl << std::endl;
        cpp << cppContent;
        cpp.close();

        std::string buildDir = m_workingDirectory + "/build";
        zz::fs::Path buildDirPath(buildDir);
        if (!buildDirPath.exist()) {
            bool dirCreated = zz::os::create_directory(buildDir);
            assert(dirCreated);
        } else {
#ifdef _MSC_VER
            zz::os::remove_dir(buildDir + "\\lib\\Release");
#else
            zz::os::remove_dir(buildDir + "/lib");
#endif
        }

        std::cout << "Check Ninja availability" << std::endl;

        std::string check_ninja = "ninja --version";
        int ret = std::system(check_ninja.c_str());
        bool ninjaAvailable = ret == EXIT_SUCCESS;

        if (ninjaAvailable) {
            std::cout << "Ninja available" << std::endl;
        } else {
            std::cout << "Ninja not available" << std::endl;
        }

        std::cout << "Configuring Cmake.." << std::endl;

        std::string buildCommand = "cmake -B" + buildDir + " -H" + m_workingDirectory;

        if (ninjaAvailable) {
            buildCommand = buildCommand + " -G Ninja";
        }

        ret = std::system(buildCommand.c_str());
        assert(ret == EXIT_SUCCESS && "The cmake configuration failed");

        std::cout << "Building.." << std::endl;

        buildCommand = "cmake --build " + buildDir + " --config Release";

        ret = std::system(buildCommand.c_str());
        assert(ret == EXIT_SUCCESS && "The compilation failed");

#ifdef _MSC_VER
        zz::fs::Path libDir(buildDir + "\\lib\\Release", true);
#else
        zz::fs::Path libDir(buildDir + "/lib");
#endif
        size_t attempts = 0;
        while (!(libDir.is_dir() && !libDir.empty()) && attempts < 1e6) {
            attempts++;
            assert(attempts != 1e6 && "The library file was not created.");

            using namespace std::chrono_literals;
            std::this_thread::sleep_for(5us);
        }

        shlibpp::SharedLibraryClassFactory<BaseClass>& shlibFactory = baseClassFactory.m_compiledEvaluableFactory;

        shlibFactory.extendSearchPath(buildDir + "\\lib\\Release");
        shlibFactory.extendSearchPath(buildDir + "/lib");

        std::cout << "Opening new library.." << std::endl;

        bool isLibOpen = shlibFactory.open((m_cleanName + "Lib").c_str(), (className + "Factory").c_str());

        assert(isLibOpen && shlibFactory.isValid() && "Unable to open compiled shared library.");

        baseClassFactory.m_compiledEvaluable = shlibFactory.create();

        assert(baseClassFactory.m_compiledEvaluable != nullptr && "The compiled instance is not valid.");

        std::cout << "Automatic generation completed!" << std::endl;
    }

    const std::vector<SqueezedMatrixRef>& evaluateGenerics() {
        for (size_t generic : m_generics) {
            m_expandedExpression[generic].buffer = m_expandedExpression[generic].partialExpression.evaluate();
        }

        return m_genericsRefs;
    }

    std::vector<levi::Registrar*> getDependencies() {
        std::vector<levi::Registrar*> deps;
        for (size_t i = 0; i < m_generics.size(); ++i) {
            std::vector<levi::Registrar*> genericDeps = m_expandedExpression[m_generics[i]].partialExpression.getDependencies();
            deps.insert(deps.end(), genericDeps.begin(), genericDeps.end());
        }
        return deps;
    }

    const std::ostringstream& getHelpersDeclaration() const {
        return m_helpersDeclarations;
    }

    const std::ostringstream& getCommonsDeclaration() const {
        return m_commonDeclarations;
    }

    const std::vector<std::ostringstream>& getFinalExpressions() const {
        return m_finalExpressions;
    }



};

#endif // AUTOGENERATEDHELPER_H