/**
 * General STREAM C++ interface.
 */


#include <cmath>
#include <iostream>
#include <H5Cpp.h>
#include <string>
#include <unistd.h>

// If "not debugging" is defined, then we're in
// debug mode and would like to active floating-point
// exceptions
#include <csignal>
#if !defined(NDEBUG) && defined(__linux__)
#   include <fenv.h>
#endif

#include <softlib/SOFTLibException.h>

#include "STREAM/stream.h"
//#include "DREAM/Init.h"
#include "DREAM/IO.hpp"
#include "DREAM/QuitException.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SFile.hpp"
#include "DREAM/Simulation.hpp"
#include "FVM/FVMException.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace std;

struct cmd_args {
    bool display_settings=false;
    bool print_adas=false;
    bool splash=true;
    string
        input_filename;
};

void display_settings(DREAM::Settings *s=nullptr) {
    if (s == nullptr)
        s = STREAM::SimulationGenerator::CreateSettings();

    cout << endl << "LIST OF STREAM SETTINGS" << endl;
    cout <<         "-----------------------" << endl;
    s->DisplaySettings();
}

void display_adas(STREAM::Simulation *sim) {
    sim->GetADAS()->PrintElements();
    cout << endl;
}

/**
 * Print the STREAMi command-line argument help.
 */
void print_help() {
    cout << "Syntax: streami INPUT [OPTIONS...]" << endl;
    cout << "   Run a startup simulation according to the specifications made" << endl;
    cout << "   in the file 'INPUT'." << endl << endl;

    cout << "OPTIONS" << endl;
    cout << "  -a           Print list of elements in ADAS database." << endl;
    cout << "  -h           Print this help." << endl;
    cout << "  -l           List all available settings in STREAM." << endl;
    cout << "  -s           Do not show the splash screen." << endl;
}

/**
 * Parse command-line arguments.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
struct cmd_args *parse_args(int argc, char *argv[]) {
    char c;

    struct cmd_args *a = new struct cmd_args;
    a->display_settings = false;

    while ((c = getopt(argc, argv, "ahls")) != -1) {
        switch (c) {
            case 'a':
                a->print_adas = true;
                break;
            case 'h':
                print_help();
                break;
            case 'l':
                display_settings();
                break;
            case 's':
                a->splash = false;
                break;
            case '?':
                cout << "Unrecognized option: " << optopt << endl;
                return nullptr;
        }
    }

    if (optind+1 < argc) {
        cout << "Too many trailing input arguments." << endl;
        return nullptr;
    } else if (optind == argc) {
        cout << "No input file specified." << endl;
        return nullptr;
    }

    a->input_filename = string(argv[optind]);
    return a;
}

void splash() {
    cout << endl;
    cout << R"(TBA: Fancy STREAM splash screen)" << endl;
    /*cout << R"( It's time to...)" << endl;
    cout << endl;
    cout << R"( * * ________   _____  ______ ___* *  ___  ___       )" << endl;
    cout << R"(  * \\   __   \/ __  \/  ____//   \ //   \/   \      )" << endl;
    cout << R"(    //  / //  / /_// /  /__ // /\  \/  / / /  /   *  )" << endl;
    cout << R"(   //  / //  /      /   __/// /_/  /  / / /  /   * * )" << endl;
    cout << R"( _//  /_//  /  /\  \   /__//  __  /  / / /  /     *  )" << endl;
    cout << R"( \\________/__/ \\__\_____//_///_/__//_//__/     *   )" << endl;
    cout << R"(               * *          * *      *         *     )" << endl;
    cout << R"(      * *       *       * *      * *             *   )" << endl;
    cout << R"(     * * *                                 ...baby...)" << endl;*/
    cout << endl;
}

/**
 * Handle floating point exceptions.
 */
void sig_fpe(int) {
    throw DREAM::FVM::FVMException("Floating-point error.");
}

/**
 * Handle a 'SIGQUIT' signal.
 */
void sig_quit(int) {
    throw DREAM::QuitException("The user requested execution to stop.");
}

/**
 * Program entry point.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
int main(int argc, char *argv[]) {
    int exit_code = 0;

    // Initialize the DREAM library
    dream_initialize();

    // Allow the user to press Ctrl+\ or Ctrl+Y to quit the simulation early
    PetscPopSignalHandler();
    std::signal(SIGQUIT, sig_quit);

#if !defined(NDEBUG) && defined(__linux__)
    std::signal(SIGFPE, sig_fpe);
#endif

    // Parse command-line arguments
    struct cmd_args *a = parse_args(argc, argv);
    if (a == nullptr)
        return -1;

    // Make sure that anything written to stdio/stderr is
    // written immediately
    cout.setf(ios_base::unitbuf);

    if (a->splash)
        splash();

    cout << "alpha version (commit " << STREAM_GIT_SHA1 << ")" << endl;

    // Except on NaN (but only in debug mode)
#if !defined(NDEBUG) && defined(__linux__)
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    DREAM::Simulation *sim = nullptr;
    try {
        DREAM::Settings *settings = STREAM::SimulationGenerator::CreateSettings();
        
        DREAM::SettingsSFile::LoadSettings(settings, a->input_filename);

        sim = STREAM::SimulationGenerator::ProcessSettings(settings);

        if (a->print_adas)
            display_adas(sim);

        sim->Run();
    } catch (DREAM::QuitException &ex) {
        DREAM::IO::PrintInfo(ex.what());
        exit_code = 0;
    } catch (DREAM::FVM::FVMException &ex) {
        DREAM::IO::PrintError(ex.what());
        exit_code = 1;
    } catch (SOFTLibException &ex) {
        DREAM::IO::PrintError(ex.what());
        exit_code = 2;
    } catch (H5::FileIException &ex) {
        DREAM::IO::PrintError(ex.getDetailMsg().c_str());
        exit_code = 3;
    }

    if (sim != nullptr) {
        try {
            sim->Save();
        } catch (H5::FileIException &ex) {
            DREAM::IO::PrintError(ex.getDetailMsg().c_str());
            exit_code = 4;
        }
    }

    // De-initialize the DREAM library
    dream_finalize();

    return exit_code;
}

