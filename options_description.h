/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file options_description.h
 * @brief OptionsDescription Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-07
 */

#ifndef __MISC_OPTIONS_DESCRIPTION_H_

#define __MISC_OPTIONS_DESCRIPTION_H_

#include <ostream>
#include <string>
#include <vector>


/**
 * @brief It is a class used to store and manipulate command line options.
 */
class OptionsDescription
{
    struct Option;
public:
    enum OptionType { kOptionBool, kOptionInt, kOptionDouble, kOptionString };

    friend std::ostream &operator <<(std::ostream &os, const OptionsDescription &desc);
    friend std::ostream &operator <<(std::ostream &os, const Option &option);

    void Parse(int &argc, char *argv[]);
    operator std::string() const;

    void AddOption(const std::string &long_name, const std::string &short_name,
            bool &bool_option, const std::string &description);

    void AddOption(const std::string &long_name, const std::string &short_name,
            int &int_option, const std::string &description);

    void AddOption(const std::string &long_name, const std::string &short_name,
            double &double_option, const std::string &description);

    void AddOption(const std::string &long_name, const std::string &short_name,
            std::string &string_option, const std::string &description);

private:
    void AddOption(const std::string &long_name, const std::string &short_name,
            void *pointer, OptionType type, const std::string &description, const std::string &default_value);

    struct Option
    {
        Option(const std::string &long_name, const std::string &short_name,
                void *pointer, OptionType type, const std::string &description, 
                const std::string &default_value);

        void Parse();
        operator std::string() const;

        std::string long_name;
        std::string short_name;
        void *pointer;
        OptionType type;
        std::string description;
        std::string default_value;
        std::string value;
    };

    std::vector<Option> options;
};

#endif

