import pandas as pd
import plotly
import plotly.graph_objects as go
import numpy as np

def plot_e_s(es_c_data):
    fig = go.Figure() 
    for k, v in es_c_data.items():
        fig.add_trace(
            go.Scatter(
                x= v[0][:,1],
                y= v[0][:,0],
                name=v[1],
                mode="markers",
                marker={'size':10+10*k}
            )
        )

    fig.update_layout(
        title="samples_events_phases",
        xaxis_title="Samples",
        yaxis_title="Events",
        font=dict(
            size=10,
            color="#7f7f7f"
        )
    )
    fig.update_layout(showlegend=True)
    return fig


df_new = pd.read_excel('PeakBinning.xlsx')
print(df_new.columns)

# Events from Peak
events = {}
for peak_id in df_new.PeakID.unique():
    events[peak_id] = df_new['Spectrum Number'][df_new.PeakID==peak_id].values 

# Co-occurence matrix
n_events = len(events)
co_mat_counts = np.zeros((n_events, n_events), dtype=float)
co_mat_samples = np.zeros((n_events, n_events), dtype=object)

for i in range(n_events):
    for j in range(n_events):
        if (events[i+1] is not None) and (events[j+1]  is not None): 
            v =  np.intersect1d(events[i+1], events[j+1])
            co_mat_counts[i,j] = len(v)
            co_mat_samples[i,j] = v

# Compute nuclei
nuclei = {'phase': [], 'samples':[],'counts':[] }
for i in range(n_events):
    for j in range(n_events):
        if (i != j) and co_mat_counts[i,i] != 0:
            if (co_mat_counts[i,i] == co_mat_counts[i,j]) and (co_mat_counts[i,i] == co_mat_counts[j,j]):
                nuclei['phase'].append((i+1,j+1))
                nuclei['samples'].append(co_mat_samples[i,j])
                nuclei['counts'].append(co_mat_counts[i,j])
                print("phase: ",(i+1,j+1), ", samples: ",co_mat_samples[i,j],", count: ",co_mat_counts[i,j])

ps = []
for val in set(nuclei['counts']):
    p = set(np.array(nuclei['phase'])[np.array(nuclei['counts'])==val].flatten())
    s = set(np.hstack(np.array(nuclei['samples'])[np.array(nuclei['counts'])==val]))
    ps.append((p,s))
    print("phase: ",p,", samples: ", s,", count: ", val)

# Plot phases in sample vs event graph
es_data = np.zeros((len(np.hstack([e for i, e in events.items() ])),3))
print(np.shape(es_data))
indx = 0 
for i, event in events.items():
    l = len(event)
    start = indx; end = indx + l
    es_data[start:end,0] = i
    es_data[start:end,1] = event
    indx += l

color_increment = round(1/float(len(ps)),2)
c = 0
es_c_data = {}
for p,s in ps:
    c = c + color_increment
    es_c_data[c] = [np.empty((0,2)),f"Phase:{p}"]
    for pp in p:
        es_c_data[c][0] = np.vstack((es_c_data[c][0],es_data[es_data[:,0]==pp ,:2]))
        es_data[es_data[:,0]==pp , 2] = c
es_c_data[0] = [es_data[es_data[:,2]==0.0 , :2] ,"No Phase"]

fig = plot_e_s(es_c_data)
fig.show()
plotly.offline.plot(fig, filename='samples_events_phases.html')