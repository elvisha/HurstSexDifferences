function [motion_seg] = motion_by_segment(scan_motion,t_start,t_end)
%compute average motion over entire segment that you're using to compute
%mean Hurst over

%scan_motion is just mxn framewise displacement where m is number of
%subjects, n is total time points (i.e. 1201)
%t_start is starting time point for segments used to calculate hurst
%t_end is last time point for segments used to calculate mean hurst

motion_seg=mean(scan_motion(:,t_start:t_end),2);

%motion_seg(:,1)=mean(scan_motion(:,11:160),2);
%motion_seg(:,2)=mean(scan_motion(:,159:308),2);
%motion_seg(:,3)=mean(scan_motion(:,307:456),2);
%motion_seg(:,4)=mean(scan_motion(:,455:604),2);
%motion_seg(:,5)=mean(scan_motion(:,603:752),2);
%motion_seg(:,6)=mean(scan_motion(:,751:900),2);
%motion_seg(:,7)=mean(scan_motion(:,899:1048),2);
%motion_seg(:,8)=mean(scan_motion(:,1047:1196),2);

end

