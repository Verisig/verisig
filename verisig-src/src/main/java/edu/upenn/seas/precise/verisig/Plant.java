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

import java.util.List;
import java.util.Map;

public class Plant {
    public List<String> states;
    public Map<Integer, Map<String, Object>> modes;
    public Map<String, Map<Integer, List<Map<String, Object>>>> glue;
    public Map<Integer, String> nameMap;

    public Plant(List<String> states, Map<Integer, Map<String, Object>> modes, Map<String, Map<Integer, List<Map<String, Object>>>> glue, Map<Integer, String> nameMap) {
        this.states = states;
        this.modes = modes;
        this.glue = glue;
        this.nameMap = nameMap;
    }
}
