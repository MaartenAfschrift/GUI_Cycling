function [casadi_path] = install_casadi()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% check casadi installation
% try to get casadi path

try casadi.GlobalOptions.getCasadiPath();
    casadi_path = casadi.GlobalOptions.getCasadiPath();
    disp(['casadi installation found in ', casadi_path]);
catch
    % first try to find casadi path in computer path
    filepath = fileparts(which('casadi-cli.exe'));
    if isempty(filepath)
        % download and install casadi if needed
        disp(' Casadi not found. Install casadi started ');

        % check operating system and matlab version
        osStruct = computer;

        % Mac operation system      
        if contains(osStruct,'MAC')
            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
            [status, cmdout] = system('uname -m');
            % Get MATLAB version
            matlabVer = version('-release');  % e.g., '2023b', '2018b'
            matlabVerNum = sscanf(matlabVer, '%d');

            if status ~= 0
                % assumping running on 	Mac classic (High Sierra or above)
                casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
            else
                machineArch = strtrim(cmdout);
                if strcmp(machineArch, 'x86_64')
                    info = 'Intel Mac';
                    casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
                elseif strcmp(machineArch, 'arm64')
                    [~, procArch] = system('arch');  % returns 'arm64' if native, 'i386' or 'x86_64' if Rosetta
                    procArch = strtrim(procArch);
                    if strcmp(procArch, 'arm64')
                        if matlabVerNum >= 2023
                            info = 'Apple Silicon (native, R2023b or later)';
                            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx_arm64-matlab2018b.zip';
                        else
                            info = 'Apple Silicon (native, older MATLAB)';
                            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
                        end
                    elseif strcmp(procArch, 'i386') || strcmp(procArch, 'x86_64')
                        if matlabVerNum >= 2018
                            info = 'Apple Silicon (Rosetta, R2018b or later)';
                            casadi_url ='https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
                        else
                            info = 'Apple Silicon (Rosetta, older MATLAB)';
                            disp('no casadi support');
                        end
                    else
                        info = ['Apple Silicon (unknown process arch: ', procArch, ')'];
                        casadi_url  = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';
                    end
                else
                    info = ['Unknown machine architecture: ', machineArch];
                    casadi_url  = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-osx64-matlab2018b.zip';

                end

            end
        % windows
        elseif contains(osStruct,'PC')
            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-windows64-matlab2018b.zip';
        % Linux
        elseif contains(osStruct,'GLN')
            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-linux64-matlab2018b.zip';
        else
            % assuming
            casadi_url = 'https://github.com/casadi/casadi/releases/download/3.7.0/casadi-3.7.0-linux64-matlab2018b.zip';
        end




        install_path = fullfile(pwd,'casadi_install');
        if ~isfolder(install_path)
            mkdir(install_path);
        else
            % all_files = dir(fullfile(install_path, '**\*.*'));
            % for i =1:length(all_files)
            %     if ~all_files(i).isdir
            %         delete(fullfile(all_files(i).folder, all_files(i).name));
            %     end
            % end
            % rmdir(install_path, 's');
            % mkdir(install_path);
        end
        websave(fullfile(install_path,'casadi_installer.zip'),casadi_url);
        disp(' unzip casadi installer');
        unzip(fullfile(install_path,'casadi_installer.zip'),install_path);
        delete(fullfile(install_path,'casadi_installer.zip'));
        disp([' Casadi installed in folder ', install_path]);
        addpath(genpath(install_path));
    else
        addpath(genpath(filepath));
    end
end