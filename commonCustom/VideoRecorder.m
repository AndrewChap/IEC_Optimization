classdef VideoRecorder < handle
    
    
% Drew Chap 2013-2017
%
% before simulation, set up recorder
%        videoRecord = true;
%        if videoRecord
%            videoRecorder = VideoRecorder(mfilename, pwd);
%            videoRecorder.zoom = '-m1';
%        end
    
% when you want to record a video frame
%        if videoRecord
%           videoWrite(videoRecorder,frameNumber);
%        end
    
% video must be closed before it can be played back
%       if videoRecord
%           closeVideo(videoRecorder);   %close video file
%       end

    properties (Access=public)
        writerObj
        imagePath
        frameRate = 30;
        videoQuality = 100;
        zoom = '-m1';
    end
    methods
        function obj = VideoRecorder(varargin)
%             keyboard
            filetitle = varargin{1};
            videoPath = varargin{2};
            fullTitle = false;
            for v = 1:length(varargin)-2
%                 if numel(varargin{v}==4)
%                     cropDimensions = varargin{v};
%                 else
                    switch varargin{v}
                        case {'-fullTitle'}
                            fullTitle = true;
                    end
%                 end
            end
            
            if fullTitle
                nameTag = filetitle;
            else
                nameTag = [filetitle,'__',clockstring];
            end
            
            mkdir('frameImages')
%             if ~exist('figure_folders','dir')
%                 mkdir('figure_folders')
%                 disp(['creating a folder ''\figure_folders'' in ',videoPath])
%             end
%             addpath([videoPath,'\figure_folders'])
%             newfolder = ['figs_',nameTag];
%             obj.imagePath = [videoPath,'\figure_folders\',newfolder,'\',nameTag];
            obj.imagePath = ['frameImages\',nameTag];
%             cd('figure_folders')
%             mkdir(newfolder)
%             cd('..')
%             if ~exist('videos','dir')   %check if the 'videos' folder in the current videoPath does not exist
%                 mkdir(videoPath,'videos')
%                 disp(['Created new folder ''videos'' in ',videoPath,'... That''s where your simulation videos will be.'])
%             end
            disp(['Creating new .avi video in ',videoPath,'\videos'])
            video_string = [nameTag];
            obj.writerObj = VideoWriter(video_string,'Motion JPEG AVI');
            set(obj.writerObj, 'FrameRate',obj.frameRate,'Quality',obj.videoQuality);
            open(obj.writerObj);
        end

        function videoWrite(obj,frameNumber)
            filename = [obj.imagePath,'_frame',num2str(frameNumber,'%05d')];  %give it a unique filename
%             tic
%             export_fig_crop_mod(filename,'-q101','-nocrop',obj.zoom,'-opengl');    %export as a png image
            export_fig_crop_mod(filename,'-q101',obj.zoom,'-opengl');    %export as a png image
%             export_fig_custom_crop(filename,'-q101',obj.zoom,'-opengl',cropDimensions);    %export as a png image

%             keyboard
            frame=imread([filename,'.png']);        % read in the image
            writeVideo(obj.writerObj,frame);        % write it to a video
%             disp(['wrote in ',num2str(toc)])
        end
        function closeVideo(obj)
            close(obj.writerObj);
        end
    end
end
        