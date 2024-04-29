%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 3 - CREATION OF TEMPLATES FOR MULTILEVEL DI ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

main_path = '/home/iesteves/FC';
template_path = [main_path, '/files/template_matrices'];
sFC_path = [main_path, '/data/results/sFC_Pearson/MIGcomplete-HCcomplete_SchaeferSubCRB7100_nonzero-50']; 
pipelinelog_path = [main_path, '/logs/pipeline_logs'];

managefolders(template_path, 'create');

pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';

% load Pearson correlation variables (iFC_all - areas x areas x subject_session, ordered by session type)
FC_file = load([sFC_path, '/Pearson-results-', atlas,'-', pipeline,'-bysubject.mat']);
nr_areas = FC_file.N_areas;
nr_subjects = FC_file.N_subjects;
filenames = FC_file.filenames_bysubject;

% log file name
DItemplates_log_filename = [pipelinelog_path,'/step3_DI_templates-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%% Generate IA and IB template matrices 
within_session_aux = zeros(nr_subjects, nr_subjects);
within_subject_aux = zeros(nr_subjects, nr_subjects);
within_group_aux = zeros(nr_subjects, nr_subjects);
within_menstrual_session_aux = zeros(nr_subjects, nr_subjects);
between_group_within_session_aux = zeros(nr_subjects, nr_subjects);

for j = 1:nr_subjects
    
    for t = 1:nr_subjects
        
        % get file 1 and split to get subject and session
        filename1 = filenames{j};
        filename1_aux = strsplit(filename1, '/');
        session_ind1 = contains(filename1_aux, 'ses-');
        subject_ind1 = contains(filename1_aux, 'sub-');
        subject1 = filename1_aux{subject_ind1};
        session1 = filename1_aux{session_ind1};
        
        % get file 2 and split to get subject and session
        filename2 = filenames{t};
        filename2_aux = strsplit(filename2, '/');
        session_ind2 = contains(filename2_aux, 'ses-');
        subject_ind2 = contains(filename2_aux, 'sub-');
        subject2 = filename2_aux{subject_ind2};
        session2 = filename2_aux{session_ind2};
        
        % same session
        if strcmp(session1, session2)
            within_session_aux(j,t) = 1;
        end
        
        % same subject
        if strcmp(subject1, subject2)
            within_subject_aux(j,t) = 1;
        end
        
        % same group
        if strcmp(subject1(1:end-3), subject2(1:end-3))
            within_group_aux(j,t) = 1;
        end
        
        % different groups, same menstrual cycle phase
        cond1= any(strcmp(session1, {'ses-preictal', 'ses-ictal', 'ses-postictal'})) & any(strcmp(session2, {'ses-premenstrual'}));
        cond2= any(strcmp(session2, {'ses-preictal', 'ses-ictal', 'ses-postictal'})) & any(strcmp(session1, {'ses-premenstrual'}));
        cond3= any(strcmp(session1, {'ses-interictal'})) & any(strcmp(session2, {'ses-midcycle'}));
        cond4= any(strcmp(session1, {'ses-midcycle'})) & any(strcmp(session2, {'ses-interictal'}));
        if cond1 || cond2 || cond3 || cond4
            between_group_within_session_aux(j,t) = 1;
        end
        
        % same menstrual cycle phase
        cond1= any(strcmp(session1, {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-premenstrual'})) & any(strcmp(session2,  {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-premenstrual'}));
        cond2= any(strcmp(session1, {'ses-midcycle', 'ses-interictal'})) & any(strcmp(session2, {'ses-midcycle', 'ses-interictal'}));
        if cond1 || cond2 || cond3 || cond4
            within_menstrual_session_aux(j,t) = 1;
        end
        
    end
    
end

%% code diagonal and upper triangular part of the template matrices as 0.5 to ignore these values
triu_ind = ones(nr_subjects, nr_subjects);
triu_ind = triu(triu_ind); 

within_group = within_group_aux;
within_group(triu_ind==1) = 0.5;

within_subject = within_subject_aux;
within_subject(triu_ind==1) = 0.5; 

within_session = within_session_aux;
within_session(triu_ind==1) = 0.5; 

within_menstrual_session = within_menstrual_session_aux;
within_menstrual_session(triu_ind==1) = 0.5; 

between_group_within_session = between_group_within_session_aux;
between_group_within_session(triu_ind==1) = 0.5; 

save([template_path, '/within_group'], 'within_group')
save([template_path, '/within_subject'], 'within_subject')
save([template_path, '/within_session'], 'within_session')
save([template_path, '/within_menstrual_session'], 'within_menstrual_session')
save([template_path, '/between_group_within_session'], 'between_group_within_session')

% Log file
DItemplates_log_msg = ['------ ', 'DI templates for: ', pipeline, '-', atlas];
disp(DItemplates_log_msg);
logCustom(DItemplates_log_msg, DItemplates_log_filename)

%% Generate IA and IB template matrices - more specific
wgroup_corr_aux = zeros(nr_subjects, nr_subjects);
wgroup_bsession_corr_aux =  zeros(nr_subjects, nr_subjects);
wsession_corr_aux =  zeros(nr_subjects, nr_subjects);
bgroup_corr_aux = zeros(nr_subjects, nr_subjects);

for j = 1:nr_subjects
    
    for t = 1:nr_subjects
        
        % get file 1 and split to get subject and session
        filename1 = filenames{j};
        filename1_aux = strsplit(filename1, '/');
        session_ind1 = contains(filename1_aux, 'ses-');
        subject_ind1 = contains(filename1_aux, 'sub-');
        subject1 = filename1_aux{subject_ind1};
        session1 = filename1_aux{session_ind1};
        
        % get file 2 and split to get subject and session
        filename2 = filenames{t};
        filename2_aux = strsplit(filename2, '/');
        session_ind2 = contains(filename2_aux, 'ses-');
        subject_ind2 = contains(filename2_aux, 'sub-');
        subject2 = filename2_aux{subject_ind2};
        session2 = filename2_aux{session_ind2};
     
        if strcmp(session1,'ses-preictal')
            if strcmp(session2, 'ses-preictal')
                wsession_corr_aux(j, t) = 1;
            elseif strcmp(session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) =  1;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) =  2;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  3;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  1;
            end
            
        elseif strcmp(session1,'ses-ictal')
            if strcmp(session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 1;
            elseif strcmp(session2, 'ses-ictal')
                wsession_corr_aux(j, t) = 2;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) =  4;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  5;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  2;
            end
        
        elseif strcmp(session1,'ses-postictal')
            if strcmp (session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 2;
            elseif strcmp (session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) = 4;
            elseif strcmp(session2, 'ses-postictal')
                wsession_corr_aux(j, t) = 3;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  6;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  3;
            end
        
        elseif strcmp(session1,'ses-interictal') 
            if strcmp (session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 3;
            elseif strcmp(session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) = 5;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) = 6;
            elseif strcmp(session2, 'ses-interictal')
                wsession_corr_aux(j, t) = 4;
            elseif strcmp(session2, 'ses-midcycle')
                bgroup_corr_aux(j, t) =  4;
            end
        
        elseif strcmp(session1,'ses-premenstrual') 
            if strcmp(session2, 'ses-premenstrual')
                wsession_corr_aux(j, t) = 5;
            elseif strcmp(session2, 'ses-midcycle')
                wgroup_bsession_corr_aux(j, t) =  7;
            elseif strcmp(session2, 'ses-preictal')
                bgroup_corr_aux(j, t) =  1;
            elseif strcmp(session2, 'ses-ictal')
                bgroup_corr_aux(j, t) =  2;
            elseif strcmp(session2, 'ses-postictal')
                bgroup_corr_aux(j, t) =  3;
            end
        
        elseif strcmp(session1,'ses-midcycle')
            if strcmp(session2, 'ses-premenstrual')
                wgroup_bsession_corr_aux(j, t) =  7;
            elseif strcmp(session2, 'ses-midcycle')
                wsession_corr_aux(j, t) = 6;
            elseif strcmp(session2, 'ses-interictal')
                bgroup_corr_aux(j, t) =  4;
            end
        end

        % same group
        if strcmp(subject1(1:end-3), subject2(1:end-3)) & strcmp(subject1(1:end-3), 'sub-patient')
            wgroup_corr_aux(j,t) = 1;
        elseif  strcmp(subject1(1:end-3), subject2(1:end-3)) & strcmp(subject1(1:end-3), 'sub-control')
            wgroup_corr_aux(j, t) = 2;
        end
        
        
    end
    
end

%% code diagonal and upper triangular part of the template matrices as 0.5
triu_ind = ones(nr_subjects, nr_subjects);
triu_ind = triu(triu_ind) ; 

wgroup_corr = wgroup_corr_aux;
wgroup_corr(triu_ind==1) = 0.5;

wgroup_bsession_corr =wgroup_bsession_corr_aux;
wgroup_bsession_corr(triu_ind==1) = 0.5 ;

wsubject_bsession_corr_aux = within_subject.*wgroup_bsession_corr;
wsubject_bsession_corr = wsubject_bsession_corr_aux;
wsubject_bsession_corr(triu_ind==1) = 0.5 ;

wsession_corr = wsession_corr_aux;
wsession_corr(triu_ind==1) = 0.5;

bgroup_corr = bgroup_corr_aux;
bgroup_corr(triu_ind==1) = 0.5;
 
save([template_path, '/wgroup_corr'], 'wgroup_corr')
save([template_path, '/wgroup_bsession_corr'], 'wgroup_bsession_corr')
save([template_path, '/wsubject_bsession_corr'], 'wsubject_bsession_corr')
save([template_path, '/wsession_corr'], 'wsession_corr')
save([template_path, '/bgroup_corr'], 'bgroup_corr')

% Log file
DItemplates_log_msg = ['------ ', 'Correlation DI templates for: ', pipeline, '-', atlas];
disp(DItemplates_log_msg);
logCustom(DItemplates_log_msg, DItemplates_log_filename)
