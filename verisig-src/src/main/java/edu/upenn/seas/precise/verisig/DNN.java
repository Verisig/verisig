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

import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class DNN {
    public Map<Integer, List<List<Double>>> weights;
    public Map<Integer, List<Double>> offsets;
    public Map<Integer, String> activations;

    public DNN(Map<Integer, List<List<Double>>> weights, Map<Integer, List<Double>> offsets, Map<Integer, String> activations) {
        this.weights = weights;
        this.offsets = offsets;
        this.activations = activations;
    }

    public int numLayers() {
        int count = 0;


        List<String> activationNames = Arrays.asList("Sigmoid", "Tanh", "Relu");
        for (String activation : activations.values()) {
            if(activationNames.contains(activation)) {
                count++;
            }

            count++;
        }

        return count;
    }

    public int numStates() {
        int count = 0;

        for( List<Double> offset : offsets.values() ) {
            if(offset.size() > count) {
                count = offset.size();
            }
        }

        return count;
    }
}
