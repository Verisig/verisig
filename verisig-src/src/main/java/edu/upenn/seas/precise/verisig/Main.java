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

import com.picklingtools.pythonesque.Tab;
import com.verivital.hyst.main.Hyst;
import edu.upenn.seas.precise.verisig.utils.SystemCommandExecutor;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.yaml.snakeyaml.Yaml;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class Main {
    public static final String VERSION = "0.9";

    @Option(name="--version", usage="print version")
    private boolean version;

    @Option(name="--help", aliases={"-h"}, usage="print help message")
    private boolean help;

    @Option(name="--verbose", aliases={"-v"}, usage="enable verbose printing")
    private boolean verbose;

    @Option(name="--spaceex-config", aliases={"-sc"}, usage="SpaceEx 'cfg' file", metaVar="PATH")
    private String spaceExConfig = null;

    @Option(name="--verisig-config", aliases={"-vc"}, usage="Verisig yaml config", metaVar="PATH")
    private String verisigConfig = null;

    @Option(name="--output-model", aliases={"-o"}, usage="Enable flow* model output")
    private boolean outputModel;

    @Option(name="--output-model-name", aliases={"-of"}, usage="Flow* model output location (overrides -o)", metaVar="PATH")
    private String outputModelName = null;

    @Option(name="--no-flowstar", aliases={"-nf"}, usage="System composition only")
    private boolean noFlowstar;

    @Option(name="--flowstar-cmd", usage="Custom flow* command", metaVar="CMD")
    private String flowstarCommand = "flowstar";

    @Argument
    private List<String> arguments = new ArrayList<>();

    public static void main(String[] args) {
        new Main().doMain(args);
    }

    public void doMain(String[] args) {
        CmdLineParser parser = new CmdLineParser(this);
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);

            if(help) {
                System.out.println("Verisig v" + VERSION + " Usage: ");
                System.out.println("verisig [options...] spaceex_xml dnn_yml");
                parser.printUsage(System.out);
                System.out.println();
                System.exit(0);
            }

            if(version) {
                System.out.println("Verisig v" + VERSION);
                System.exit(0);
            }

            if( arguments.size() < 2 ) {
                throw new CmdLineException(parser, "Not enough arguments given");
            }
        } catch( CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("verisig [options...] spaceex_xml dnn_yml");
            parser.printUsage(System.err);
            System.err.println();

            return;
        }

        Path inputModelPath = Paths.get(arguments.get(0));
        Path dnnPath = Paths.get(arguments.get(1));

        if(!inputModelPath.toString().endsWith(".xml")) {
            System.err.println("Input model should be SpaceEx XML file: " + inputModelPath);
            return;
        }

        if(spaceExConfig == null) {
            spaceExConfig = inputModelPath.toString().replace(".xml", ".cfg");
            if(verbose) System.out.println("Assuming SpaceEx config of: " + spaceExConfig);
        }

        String[] hystArgs;
        if( verbose ) {
            System.out.println("Converting SpaceEx model using HyST");
            hystArgs = new String[]{"-v", "-tool", "verisig", "-in-memory", "-passes", "complete_flow", "\"\"", "-input", inputModelPath.toString(), spaceExConfig};
        } else {
            hystArgs = new String[]{"-tool", "verisig", "-in-memory", "-passes", "complete_flow", "\"\"", "-input", inputModelPath.toString(), spaceExConfig};
        }

        int returnCode = Hyst.runWithArguments(hystArgs);
        if( returnCode != Hyst.ExitCode.SUCCESS.ordinal() ) {
            System.err.println("Hyst conversion was unsuccessful");
            System.exit(returnCode);
        }

        if( verbose ) {
            System.out.println("SpaceEx model successfully converted");
        }

        Plant plant = loadPlant();
        DNN dnn;

        if(verbose) {
            System.out.println("Loading DNN from:" + dnnPath);
        }

        try{
            dnn = loadDNN(dnnPath.toString());
        } catch (Exception e) {
            System.err.println("Error loading dnn from: " + dnnPath.toString());
            e.printStackTrace();
            return;
        }

        if(verbose) {
            System.out.println("DNN successfully loaded");
        }

        if( verisigConfig == null ) {
            verisigConfig = inputModelPath.toString().replace(".xml", ".yml");
            if(verbose) System.out.println("Assuming Verisig config of: " + verisigConfig);
        }


        if(verbose) {
            System.out.println("Loading Verisig configuration file from: " + verisigConfig);
        }
        Map<String, Object> config;
        try {
            config = loadConfig(verisigConfig);
        } catch (Exception e) {
            System.err.println("Error loading config from: " + verisigConfig);
            e.printStackTrace();
            return;
        }
        if(verbose) {
            System.out.println("Verisig configuration file successfully loaded");
        }

        ComposedHybridSystem system = new ComposedHybridSystem(config, dnn, plant);

        try {
            ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();

            if(verbose) {
                System.out.println("Composing DNN and hybrid system");
            }
            system.write(new BufferedWriter(new OutputStreamWriter(byteArrayOutputStream)));
            if(verbose) {
                System.out.println("Composed system successful");
            }

            if( outputModel || outputModelName != null) {
                if( outputModelName == null ) {
                    outputModelName = inputModelPath.toString().replace(".xml", ".model");
                    if(verbose) System.out.println("Assuming flow* model output of: " + outputModelName);
                }

                FileWriter fileWriter = new FileWriter(outputModelName);
                fileWriter.write(byteArrayOutputStream.toString());
                fileWriter.close();
            }

            if( !noFlowstar ) {
                if(verbose) {
                    System.out.println("Running flowstar command: " + flowstarCommand + " on composed system");
                }

                List<String> command = new ArrayList<>();
                command.add(flowstarCommand);

                if( system.plottingEnabled ) {
                    command.add("-p");
                }

                if( system.dumpingEnabled ) {
                    command.add("-d");
                }

                command.add(dnnPath.toString());
                SystemCommandExecutor executor = new SystemCommandExecutor(command, true);
                executor.executeCommand(byteArrayOutputStream.toString());

                if(verbose) {
                    System.out.println("Flowstar command finished");
                }
            }

        } catch(Exception e) {
            e.printStackTrace();
            //
        }
    }

    protected static Plant loadPlant() {
        VerisigPrinter printer = VerisigPrinter.getInstance();
        Tab plant = printer.plant;

        if( !plant.containsKey("states") || !plant.containsKey("modes")
                || !plant.containsKey("glue") || !plant.containsKey("name_map") ) {
            throw new RuntimeException("Plant is not formatted correctly.");
        }

        return new Plant((List<String>)plant.get("states"),
                (Map<Integer, java.util.Map<String, Object>>)plant.get("modes"),
                (Map<String, Map<Integer, List<Map<String, Object>>>>)plant.get("glue"),
                (Map<Integer, String>)plant.get("name_map"));
    }

    protected static DNN loadDNN(String filename) throws IOException {
        Yaml yaml = new Yaml();
        try( FileInputStream stream = new FileInputStream(filename)) {
            Map<String, Object> yml = yaml.load(stream);

            if( !yml.containsKey("weights") || !yml.containsKey("offsets") || !yml.containsKey("activations")) {
                throw new RuntimeException("DNN is not formatted correctly.");
            }
            return new DNN((Map<Integer, List<List<Double>>>)yml.get("weights"),
                    (Map<Integer, List<Double>>) yml.get("offsets"),
                    (Map<Integer, String>) yml.get("activations"));
        }
    }

    protected static Map<String, Object> loadConfig(String filename) throws IOException {
        Yaml yaml = new Yaml();
        try( FileInputStream stream = new FileInputStream(filename)) {
            return yaml.load(stream);
        }
    }
}
