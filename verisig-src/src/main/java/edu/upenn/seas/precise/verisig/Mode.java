/*
Copyright (C) 2019 Radoslav Ivanov, Taylor J Carpenter, James Weimer, Rajeev Alur, George J. Pappa, Insup Lee

This file is part of Verisig.

Verisig is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Verisig is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Verisig.  If not, see <https://www.gnu.org/licenses/>.
*/
package edu.upenn.seas.precise.verisig;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Mode {
    public String name;
    public String odetype;
    public Map<String, String> dynamics;
    public List<String> invariant;

    public Mode(String name) {
        this(name, new HashMap<>(), new ArrayList<>());
    }

    public Mode(String name, List<String> invariant) {
        this(name, new HashMap<>(), invariant);
    }

    public Mode(String name, Map<String, String> dynamics) {
        this(name, dynamics, new ArrayList<>());
    }

    public Mode(String name, Map<String, String> dynamics, List<String> invariant) {
        this.name = name;
        this.odetype = "nonpoly ode";
        this.dynamics = dynamics;
        this.invariant = invariant;
    }


}
