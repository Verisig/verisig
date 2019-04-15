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
import java.util.List;

public class Jump {
    public String fromMode;
    public String toMode;
    public List<String> guard;
    public List<String> reset;

    public Jump(String fromMode, String toMode) {
        this(fromMode, toMode, new ArrayList<>(), new ArrayList<>());
    }

    public Jump(String fromMode, String toMode, List<String> guard, List<String> reset) {
        this.fromMode = fromMode;
        this.toMode = toMode;
        this.guard = guard;
        this.reset = reset;
    }
}
