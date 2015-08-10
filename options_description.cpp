/**
 * @file options_description.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-07
 */

#include "options_description.h"

#include <getopt.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

ostream &operator <<(ostream &os, const OptionsDescription &desc) {
    for (unsigned i = 0; i < desc.options.size(); ++i)
        os << desc.options[i] << endl;

    return os;
}

ostream &operator <<(ostream &os, const OptionsDescription::Option &option) {
    return os << string(option);
}

void OptionsDescription::Parse(int &argc, char *argv[]) {
    struct option long_options[options.size() + 1];
    string short_options;

    for (unsigned i = 0; i < options.size(); ++i) {
        long_options[i].name = options[i].long_name.c_str();
        long_options[i].flag = 0;
        long_options[i].val =
            (options[i].short_name != "" ? options[i].short_name[0] : 0);

        if (options[i].type == kOptionBool) {
            if (options[i].short_name != "")
                short_options += options[i].short_name;

            long_options[i].has_arg = no_argument;
        }
        else {
            if (options[i].short_name != "")
                short_options += options[i].short_name + ":";

            long_options[i].has_arg = required_argument;
        }
    }

    long_options[options.size()].name = 0;
    long_options[options.size()].has_arg = 0;
    long_options[options.size()].flag = 0;
    long_options[options.size()].val = 0;

    while (true) {
        int index = -1;
        int ch = getopt_long(argc, argv, short_options.c_str(), long_options, &index);

        if (ch == -1)
            break;

        if (ch == '?')
            throw logic_error("uknown option");

        //throw exception();

        if (ch != 0) {
            string s;
            s += ch;

            for (unsigned i = 0; i < options.size(); ++i) {
                if (options[i].short_name == s) {
                    index = i;
                    break;
                }
            }
        }

        if (options[index].type == kOptionBool)
            options[index].value = "true";
        else
            options[index].value = optarg;
    }

    int index = 1;

    for (int i = optind; i < argc; ++i)
        argv[index++] = argv[i];

    argc = index;

    for (unsigned i = 0; i < options.size(); ++i)
        options[i].Parse();
}

OptionsDescription::operator string() const {
    stringstream ss;

    for (unsigned i = 0; i < options.size(); ++i)
        ss << options[i] << endl;

    return ss.str();
}

void OptionsDescription::AddOption(const std::string &long_name, const std::string &short_name,
                                   bool &bool_option, const std::string &description) {
    AddOption(long_name, short_name, &bool_option, kOptionBool, description, "");
}

void OptionsDescription::AddOption(const std::string &long_name, const std::string &short_name,
                                   int &int_option, const std::string &description) {
    stringstream ss;
    ss << int_option;
    AddOption(long_name, short_name, &int_option, kOptionInt, description, ss.str());
}

void OptionsDescription::AddOption(const std::string &long_name, const std::string &short_name,
                                   double &double_option, const std::string &description) {
    stringstream ss;
    ss << double_option;
    AddOption(long_name, short_name, &double_option, kOptionDouble, description, ss.str());
}

void OptionsDescription::AddOption(const std::string &long_name, const std::string &short_name,
                                   std::string &string_option, const std::string &description) {
    AddOption(long_name, short_name, &string_option, kOptionString, description, string_option);
}

void OptionsDescription::AddOption(const string &long_name, const string &short_name,
                                   void *pointer, OptionType type, const string &description, const string &default_value) {
    for (unsigned i = 0; i < options.size(); ++i) {
        if (options[i].long_name == long_name
                || (short_name != "" && options[i].short_name == short_name)) {
            throw logic_error("option already exist: " + long_name + ", " + short_name);
        }
    }

    options.push_back(Option(long_name, short_name, pointer, type, description, default_value));
}

OptionsDescription::Option::Option(const string &long_name, const string &short_name,
                                   void *pointer, OptionType type, const string &description,
                                   const string &default_value) {
    if (short_name.size() > 1 || long_name.size() < 2) {
        throw logic_error("invalid option: " + long_name + ", " + short_name);
    }

    this->long_name = long_name;
    this->short_name = short_name;
    this->pointer = pointer;
    this->type = type;
    this->description = description;
    this->default_value = default_value;
}

void OptionsDescription::Option::Parse() {
    if (value == "")
        value = default_value;

    if (value != "") {
        switch (type) {
        case kOptionBool:
            if (value != "" || value == "false")
                *(bool *)pointer = true;

            break;

        case kOptionInt:
            *(int *)pointer = atoi(value.c_str());
            break;

        case kOptionDouble:
            *(double *)pointer = atof(value.c_str());
            break;

        case kOptionString:
            *(string *)pointer = value;
        }
    }
}

OptionsDescription::Option::operator string() const {
    string s;

    if (short_name != "")
        s += "  -" + short_name + ", ";
    else
        s += "      ";

    s += "--" + long_name;

    if (type != kOptionBool)
        s += " arg";

    if (default_value != "")
        s += " (=" + default_value + ")";

    while (s.size() < 40)
        s += " ";

    s += " " + description;

    return s;
}

