""" Prepare Bloch Sphere settings """
db=Bloch()
colors = 'b'
db.point_color = colors
db.point_marker = ['o']

# qsave(result_states, "anim/odd")
result_states = qload("anim/odd")

# Animation by saving looooots of figs 
for i in range(0,len(result_states), 10):
    db.add_states(result_states[i].ptrace(0), kind='point')
    db.render()
    db.fig.savefig("anim/odd_circle_%03d.png" % (i/10))

db.add_states(result_states[-1].ptrace(0))
db.render()
db.fig.savefig("anim/odd_circle_final.png")

# ffmpeg -framerate 50 -i odd_circle_%3d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

