// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROBOPTIM_CORE_IPOPT_PARAMETERS_UPDATER_HH
# define ROBOPTIM_CORE_IPOPT_PARAMETERS_UPDATER_HH

# include <string>

# include <boost/variant/static_visitor.hpp>

# include <coin/IpSmartPtr.hpp>

# include <roboptim/core/function.hh>

namespace roboptim
{
  namespace
  {
    struct IpoptParametersUpdater : public boost::static_visitor<bool>
    {
      explicit IpoptParametersUpdater
      (const Ipopt::SmartPtr<Ipopt::IpoptApplication>& app,
       const std::string& key)
	: boost::static_visitor<bool> (),
	  app (app),
	  key (key)
      {}

      bool operator () (const Function::value_type& val) const
      {
	return app->Options ()->SetNumericValue (key, val);
      }

      bool operator () (const int& val) const
      {
	return app->Options ()->SetIntegerValue (key, val);
      }

      bool operator () (const std::string& val) const
      {
	return app->Options ()->SetStringValue (key, val);
      }

      bool operator () (const bool& val) const
      {
	return app->Options ()->SetStringValue (key, val? "yes" : "no");
      }

      template <typename T>
      bool operator () (const T&) const
      {
        throw std::runtime_error ("option type not supported by Ipopt.");
      }

    private:
      const Ipopt::SmartPtr<Ipopt::IpoptApplication>& app;
      const std::string& key;
    };
  } // end of anonymous namespace.
} // end of namespace roboptim.

#endif //! ROBOPTIM_CORE_IPOPT_PARAMETERS_UPDATER_HH
