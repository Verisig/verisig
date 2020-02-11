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

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ComposedHybridSystem {
    public Map<String,Object> config;
    public DNN dnn;
    public Plant plant;
    public List<String> allStates;
    public List<String> dnnStates;
    public List<String> plantStates;
    public List<Mode> modes;
    public List<Jump> jumps;
    public Map<String,String> properties;

    public boolean plottingEnabled = false;
    public boolean dumpingEnabled = false;

    public static List<String> PROPERTY_ORDER = Arrays.asList("adaptive steps", "time", "remainder estimation", "identity precondition",
            "gnuplot", "fixed orders", "cutoff", "precision", "output", "max jumps", "print");
    public static Pattern ATAN_PATTERN = Pattern.compile("_arc_(_?[^_\\W]+)_(_?[^_\\W]+)(?:_([^_\\W]+))?");
    public static Pattern SQRT_PATTERN = Pattern.compile("_sqrt_(_?[^_\\W]+)_(_?[^_\\W]+)(?:_([^_\\W]+))?");
    public static Pattern DIV_PATTERN = Pattern.compile("_div_(_?[^_\\W]+)_(_?[^_\\W]+)(?:_([^_\\W]+))?");
    public static Pattern COS_PATTERN = Pattern.compile("_cos_(_?[^_\\W]+)_(_?[^_\\W]+)(?:_([^_\\W]+))?");
    public static Pattern SIN_PATTERN = Pattern.compile("_sin_(_?[^_\\W]+)_(_?[^_\\W]+)(?:_([^_\\W]+))?");

    public ComposedHybridSystem(Map<String, Object> config, DNN dnn, Plant plant) {
        this.config = config;
        this.dnn = dnn;
        this.plant = plant;

        composeVariables();
        composeModes();
        composeJumps();

        populateProperties();

        if( config.containsKey("plot") ) {
            plottingEnabled = (boolean) config.get("plot");
        }

        if( config.containsKey("dump") ) {
            dumpingEnabled = (boolean) config.get("dump");
        }
    }

    protected void composeVariables() {
        dnnStates = new ArrayList<>();
        plantStates = new ArrayList<>();
        for( String state : plant.states ) {
            if( state.startsWith("_f") )
            {
                dnnStates.add(state);
            } else if( state.equals("clock") ) {
                // ignore
            } else {
                plantStates.add(state);
            }
        }

        allStates = new ArrayList<>(dnnStates);
        allStates.addAll(plantStates);
        allStates.add("clock");
    }

    protected void composeModes() {
        modes = createDNNModes();
        modes.addAll(createPlantModes());
    }

    protected List<Mode> createDNNModes() {
        List<Mode> modes = new ArrayList<>();
        modes.add( createDNNMode(0));
        modes.add( createDNNMode(1, "DNN"));

        return modes;
    }

    protected Mode createDNNMode(int modeIndex) {
        return createDNNMode(modeIndex, "");
    }

    protected Mode createDNNMode(int modeIndex, String prefix) {
        Mode mode = new Mode(prefix + "m" + modeIndex);

        for(String state : allStates) {
            mode.dynamics.put(state, state + "\' = 0");
        }

        mode.dynamics.put("clock", "clock\' = 1");
        mode.invariant.add("clock <= 0");

        return mode;
    }

    protected List<Mode> createPlantModes() {
        List<Mode> modes = new ArrayList<>();

        for (Map.Entry<Integer, Map<String, Object>> entry : plant.modes.entrySet()) {
            modes.add(createPlantMode(entry.getKey(), entry.getValue()));
        }

        return modes;
    }

    protected Mode createPlantMode(int modeId, Map<String,Object> modeMap) {
        String modeName = (String) modeMap.get("name");
        Map<String, String> dynamics = (Map<String, String>) modeMap.get("dynamics");
        List<String> invariant = (List<String>) modeMap.get("invariant");

        Matcher matcher = ATAN_PATTERN.matcher(modeName);
        String prefix = "_arc";

        if(!matcher.matches()) {
            matcher = SQRT_PATTERN.matcher(modeName);
            prefix = "_sqrt";

            if(!matcher.matches()) {
                matcher = DIV_PATTERN.matcher(modeName);
                prefix = "_div";

                if(!matcher.matches()) {
                    matcher = COS_PATTERN.matcher(modeName);
                    prefix = "_cos";

                    if(!matcher.matches()) {
                        matcher = SIN_PATTERN.matcher(modeName);
                        prefix = "_sin";
                    }
                }
            }
        }

        if(matcher.matches()) {
            String storeVar = matcher.group(1);
            String inputVar = matcher.group(2);


            int storeVarIndex = -1;
            int inputVarIndex = -1;

            for( int i = 0; i < allStates.size(); i++) {
                if(storeVar.equals(allStates.get(i))) {
                    storeVarIndex = i;
                }

                if(inputVar.equals(allStates.get(i))) {
                    inputVarIndex = i;
                }
            }

            if( storeVarIndex < 0 || inputVarIndex < 0) {
                throw new RuntimeException("Either storage variable or input variable of " + modeName + " is not known");
            }

            if(matcher.groupCount() == 3) {
                modeName = prefix + "_" + storeVarIndex + "_" + inputVarIndex + "_" + matcher.group(3);
            } else {
                modeName = prefix + "_" + storeVarIndex + "_" + inputVarIndex + "_1";
            }

            plant.nameMap.put(modeId, modeName);
        }

        return new Mode(modeName, dynamics, invariant);
    }

    protected void composeJumps() {
        jumps = createDNNJumps();
        jumps.addAll(createPlantJumps());
        jumps.addAll(createDNN2PlantJumps());
        jumps.addAll(createPlant2DNNJumps());
    }

    protected List<Jump> createDNNJumps() {
        Jump jump = new Jump("m0", "DNNm1");

        /*for( int state = 0; state < dnn.numStates(); state++ ) {
            jump.reset.add("_f" + (state+1) +"\' := _f" + (state+1));
        }*/

        jump.reset.add("clock\' := 0");
        List<Jump> jumps = new ArrayList<>();
        jumps.add(jump);

        return jumps;
    }

    protected List<Jump> createPlantJumps() {
        List<Jump> jumps = new ArrayList<>();

        for (Map.Entry<Integer, Map<String, Object>> entry : plant.modes.entrySet()) {
            String fromMode = plant.nameMap.get(entry.getKey());
            for (Map.Entry<Integer, List<Map<String, Object>>> transitionEntry : ((Map<Integer, List<Map<String, Object>>>) entry.getValue().get("transitions")).entrySet()) {
                String toMode = plant.nameMap.get(transitionEntry.getKey());
                for (Map<String, Object> transition : transitionEntry.getValue()) {
                    jumps.add(createPlantJump(fromMode, toMode, transition));
                }
            }
        }

        return jumps;
    }

    protected Jump createPlantJump(String fromMode, String toMode, Map<String, Object> transition) {
        return new Jump(fromMode, toMode, (List<String>)transition.get("guard"), (List<String>)transition.get("reset"));
    }

    protected List<Jump> createDNN2PlantJumps() {
        List<Jump> jumps = new ArrayList<>();

        String dnnMode = "DNNm1";
        for (Map.Entry<Integer, List<Map<String, Object>>> entry : plant.glue.get("dnn2plant").entrySet()) {
            String toMode = plant.nameMap.get(entry.getKey());

            for (Map<String, Object> transition : entry.getValue()) {
                jumps.add(createPlantJump(dnnMode, toMode, transition));
            }
        }

        return jumps;
    }

    protected List<Jump> createPlant2DNNJumps() {
        List<Jump> jumps = new ArrayList<>();

        String dnnMode = "m0";

        for (Map.Entry<Integer, List<Map<String, Object>>> entry : plant.glue.get("plant2dnn").entrySet()) {
            String fromMode = plant.nameMap.get(entry.getKey());

            for (Map<String, Object> transition : entry.getValue()) {
                jumps.add(createPlantJump(fromMode, dnnMode, transition));
            }
        }

        return jumps;
    }

    protected void populateProperties() {
        properties = new HashMap<>();

        properties.put("adaptive steps", "{min 1e-6, max 0.1}");
        properties.put("time", "100");
        properties.put("remainder estimation", "1e-1");
        properties.put("identity precondition", "");
        properties.put("gnuplot", "octagon clock, _f1");
        properties.put("fixed orders", "4");
        properties.put("cutoff", "1e-18");
        properties.put("precision", "100");
        properties.put("output", "autosig");
        properties.put("max jumps", "100");
        properties.put("print", "on");

        // Overwrite defaults with config file
        for(String key : properties.keySet()) {
            if( config.containsKey(key) ) {
                properties.put(key, String.valueOf(config.get(key)));
            }
        }
    }

    protected void writeInitialCondition(BufferedWriter stream) throws IOException {
        stream.write("\tinit\n");
        stream.write("\t{\n");

        Map<String, Object> initialization = (Map<String,Object>)config.get("init");
        stream.write("\t\t" + ((String)initialization.get("mode")).trim() + "\n");
        stream.write("\t\t{\n");

        for( String state : (List<String>)initialization.get("states") ) {
            stream.write("\t\t\t" + state.trim() + "\n");
        }

        for(String state : allStates) {
            boolean notInit = true;

            Pattern pat = Pattern.compile("([^\\d\\w_]+x[^\\d\\w_]+)|(\\sx[^\\d\\w_]+)|^x[^\\d\\w_]+".replace("x", state));
            for(String prop : (List<String>)initialization.get("states")) {
                if( pat.matcher(prop).lookingAt() ) {
                    notInit = false;
                    break;
                }
            }

            if( notInit ) {
                stream.write("\t\t\t" + state + " in [0, 0]\n");
            }
        }

        stream.write("\t\t}\n");
        stream.write("\t}\n");
    }

    protected void writeSafetyProperties(BufferedWriter stream) throws IOException {
        stream.write("unsafe\n");
        stream.write("{\n");

        List<Map<String, Object>> unsafe = (List<Map<String,Object>>) config.get("unsafe");
        for(Map<String,Object> unsafeState : unsafe) {
            stream.write("\t" + unsafeState.get("mode") + "\n");
            stream.write("\t{\n");

            for(String prop : (List<String>)unsafeState.get("states")) {
                stream.write("\t\t" + prop + "\n");
            }

            stream.write("\t}\n");
        }

        stream.write("}");
    }

    protected void write(BufferedWriter stream) throws IOException {
        stream.write("hybrid reachability\n");
        stream.write("{\n");

        stream.write("\t state var ");
        for( int i = 0; i < allStates.size(); i++ ) {
            stream.write(allStates.get(i));
            if( i < allStates.size() - 1 ) {
                stream.write(", ");
            }
        }

        stream.write("\n\n");

        stream.write("\tsetting\n");
        stream.write("\t{\n");

        for( String key : PROPERTY_ORDER ) {
            stream.write("\t\t" + key + " " + String.valueOf(properties.get(key)).trim() + "\n");
        }

        stream.write("\t}\n\n");

        stream.write("\tmodes\n");
        stream.write("\t{\n");

        for( Mode mode : modes ) {
            stream.write("\t\t" + mode.name + "\n");
            stream.write("\t\t{\n");
            stream.write("\t\t\t" + mode.odetype + "\n");
            stream.write("\t\t\t{\n");

            for( String dynamics : mode.dynamics.values() ) {
                stream.write("\t\t\t\t" + dynamics + "\n");
            }

            stream.write("\t\t\t}\n");
            stream.write("\t\t\tinv\n");
            stream.write("\t\t\t{\n");

            for( String invariant : mode.invariant ) {
                stream.write("\t\t\t\t" + invariant + "\n");
            }

            stream.write("\t\t\t}\n");
            stream.write("\t\t}\n");
        }

        stream.write("\t}\n");

        stream.write("\tjumps\n");
        stream.write("\t{\n");

        for( Jump jump : jumps ) {
            stream.write("\t\t" + jump.fromMode + " -> " + jump.toMode + "\n");
            stream.write("\t\tguard { ");

            for( String guard : jump.guard ) {
                stream.write(guard + " ");
            }

            stream.write("}\n");
            stream.write("\t\treset { ");

            for( String reset : jump.reset ) {
                stream.write(reset + " ");
            }

            stream.write("}\n");
            stream.write("\t\tinterval aggregation\n");
        }

        stream.write("\t}\n");

        writeInitialCondition(stream);

        stream.write("}\n");

        writeSafetyProperties(stream);
        stream.flush();
    }

}
