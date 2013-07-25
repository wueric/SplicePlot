import drawing
import math
import plot_settings
import argparse
import pandas

class Axis():

    def __init__(self,start_radius,end_radius,angle,true_minimum,true_maximum,step):
        self.start_radius = start_radius
        self.end_radius = end_radius
        self.angle = angle
        self.step = step
        self.minimum = true_minimum
        self.maximum = true_maximum
    
    @classmethod
    def auto_boundary_axis(cls,start_radius,end_radius,angle,minimum,maximum,step):
        low = math.floor(minimum / step) * step
        high = math.ceil(maximum / step) * step
        return cls(start_radius,end_radius,angle,low,high,step)

    def radius_map(self,value):
        assert value >= self.minimum and value <= self.maximum, "Invalid value"

        return self.start_radius + (value - self.minimum) * (self.end_radius - self.start_radius) / (self.maximum - self.minimum)

def draw_hive_plot(file_name,
                data,
                hive_plot_settings,
                genotype_ordering):

    """ Writes a .svg file containing the plot using the settings passed in as parameters

    file_name is the name of the resulting .svg
    data is a pandas.DataFrame object containing the splicing expression data
    dimension is the dimension of the resulting .svg file
    axis_angles is a list containing the angles of the axes of the plot in degrees
    axis_thickness is a number describing the thickness of the axes of the plot
    axis_color is a 3-element list containing the RGB values for the color of the axes
    axis_start_radius is the radius of the beginning of the axes, measured from the center of the plot
    axis_end_radius is the radius of the end of the axes, measured from the center of the plot
    axis_step is the step size of the axis. This controls the spacing of the tick marks on the axis
    custom_scale is a list containing lists of the lower and upper bounds for each axis
    axis_label_size is a number describing the size of the axis labels
    axis_label_radius is a list containing the radii of each of the axis labels, measured from the center of the plot
    use_ticks is a boolean which determines whether tick marks are drawn on the axes
    tick_length is a number describing the length of the tick marks. Used only if use_ticks is true
    tick_label_from_axis is a number determining the distance of the tick labels from the axis
    tick_label_size is a number determining the font size of the tick label. It is set to 0 if no tick label is drawn
    bezier_colors is a list containing lists which hold the RGB values of each of the bezier curves
    bezier_thickness is a number determining the thickness of the bezier curves
    include_key is a boolean determining whether the key is drawn
    key_title_size is a number determining the font size of the key title
    key_text_color is a 3-element list containing the RGB values for the color of the key text
    key_position is a 2-element list describing the cartesian coordinates of the key. (0,0) is the center of the plot
    key_font_size is a number determining the font size of the labels in the key
    """

    # parse hive_plot_settings
    dimension=hive_plot_settings['dimension']
    axis_angles=hive_plot_settings['axis_angles']
    axis_thickness=hive_plot_settings['axis_thickness']
    axis_color=hive_plot_settings['axis_colors']
    axis_start_radius=hive_plot_settings['axis_start_radius']
    axis_end_radius=hive_plot_settings['axis_end_radius']
    axis_step=hive_plot_settings['axis_subdivision']
    custom_scale = hive_plot_settings['custom_scale']
    axis_label_size=hive_plot_settings['axis_label_size']
    axis_label_radius=hive_plot_settings['axis_label_radius']
    use_ticks=hive_plot_settings['tick_marks']
    tick_length=hive_plot_settings['tick_height']
    tick_label_from_axis=hive_plot_settings['tick_label_distance']
    tick_label_size=hive_plot_settings['tick_label_font_size']
    bezier_colors=hive_plot_settings['bezier_colors']
    bezier_thickness=hive_plot_settings['bezier_thickness']
    include_key=hive_plot_settings['include_key']
    key_title_size=hive_plot_settings['key_title_size']
    key_text_color=hive_plot_settings['key_text_color']
    key_position=hive_plot_settings['key_position']
    key_font_size=hive_plot_settings['key_font_size']


    f1 = open(file_name,'w+')

    # set up svg
    f1.write('<?xml version="1.0"?>\n')
    f1.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')

    f1.write('<svg width="{0}in" height="{0}in" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {0} {0}">\n'.format(int(dimension)))

    f1.write('<g transform="translate({0},{0}) scale(1, -1)">\n'.format((dimension / 2)))

    # get dimensions of data frame
    data_shape = data.shape

    #######################
    # create axis objects #
    #######################
    number_of_axes = data.shape[1]-1
    axis_angles = axis_angles[0:number_of_axes]

    # determine minimum and maximum for each axis
    axis_minima = list(data.min(0)[1:])
    axis_maxima = list(data.max(0)[1:])

    all_axes = [None] * number_of_axes

    for i in range(number_of_axes):
        if custom_scale:
            all_axes[i] = Axis(axis_start_radius, axis_end_radius, axis_angles[i], custom_scale[i][0],custom_scale[i][1],axis_step)
        else:
            all_axes[i] = Axis.auto_boundary_axis(axis_start_radius,axis_end_radius,axis_angles[i],axis_minima[i],axis_maxima[i],axis_step)

    # assign colors to genotypes
    unique_genotypes = genotype_ordering
    color_assignment = {}
    for i in range(0,len(unique_genotypes)):
        color_assignment[unique_genotypes[i]] = bezier_colors[i]
        
    # draw in the bezier curves
    for individual_name in data.index:
        individual = data.loc[individual_name]

        indiv_genotype = individual[0]
        indiv_expression = individual[1:]

        for i in range(0,len(indiv_expression)-1):
            start_axis = all_axes[i]
            end_axis = all_axes[i+1]

            bezier_origin_x, bezier_origin_y = polar_to_cartesian(start_axis.radius_map(indiv_expression[i]),start_axis.angle)
            bezier_dest_x, bezier_dest_y = polar_to_cartesian(end_axis.radius_map(indiv_expression[i+1]),end_axis.angle)

            start_to_end_x, start_to_end_y = bezier_dest_x - bezier_origin_x, bezier_dest_y - bezier_origin_y

            mid_x, mid_y = 0.5 * (bezier_origin_x + bezier_dest_x), 0.5 * (bezier_origin_y + bezier_dest_y)
            control_x, control_y = 0.6 * start_to_end_y + mid_x, -0.6 * start_to_end_x + mid_y

            drawing.draw_bezier(f1,bezier_origin_x,bezier_origin_y,bezier_dest_x,bezier_dest_y,control_x,control_y,color=color_assignment[indiv_genotype],thickness=bezier_thickness)

    # draw axes and tick marks
    for axis in all_axes:
        start_x, start_y = polar_to_cartesian(axis.start_radius,axis.angle)
        end_x, end_y = polar_to_cartesian(axis.end_radius,axis.angle)

        drawing.draw_line_polar(f1,start_x,start_y,axis.end_radius-axis.start_radius,axis.angle,thickness=axis_thickness,color=axis_color)

        if use_ticks:
            i = axis.minimum
            while i <= axis.maximum:
                radius = axis.radius_map(i)
                tick_start_x, tick_start_y = polar_to_cartesian(radius,axis.angle)
                drawing.draw_line_polar(f1,tick_start_x,tick_start_y,tick_length / 2,axis.angle+90,thickness=0.3*axis_thickness,color=axis_color)
                drawing.draw_line_polar(f1,tick_start_x,tick_start_y,tick_length / 2,axis.angle-90,thickness=0.3*axis_thickness,color=axis_color)

                if tick_label_size:
                    normal_vector_x, normal_vector_y = polar_to_cartesian(tick_label_from_axis,axis.angle+90)
                    label_x, label_y = tick_start_x + normal_vector_x, normal_vector_y + tick_start_y

                    # draw right side up
                    if (axis.angle % 360 <= 90) or (axis.angle % 360 >= 270):
                        drawing.draw_text(f1,i,label_x,label_y,tick_label_size,angle=-axis.angle,color=axis_color)

                    # draw upside down
                    else:
                        drawing.draw_text(f1,i,label_x,label_y,tick_label_size,angle=-axis.angle - 180,color=axis_color)

                i += axis.step

    # draw axis labels
    for i in range(len(data.columns[1:])):
        x, y = polar_to_cartesian(axis_label_radius[i], axis_angles[i])
        drawing.draw_multiline_text(f1,data.columns[i+1].replace(':','\n'),x,y,axis_label_size,axis_color)

    # draw key
    if include_key:
        drawing.draw_text_left(f1,data.columns[0],key_position[0],key_position[1],key_title_size,color=key_text_color)

        new_reference_y = key_position[1] - key_title_size - key_font_size / 2
        spacing = key_font_size * 2

        for i in range(len(unique_genotypes)):
            genotype = unique_genotypes[i]
            color = color_assignment[genotype]

            drawing.draw_rectangle(f1,key_position[0],new_reference_y - i * spacing,key_font_size,key_font_size,color)
            drawing.draw_text_left(f1,genotype,key_position[0]+spacing,new_reference_y - i * spacing,key_font_size,color=key_text_color)
                
    # close svg
    f1.write('</g>\n')
    f1.write('</svg>\n')
    f1.close()    

def polar_to_cartesian(radius,degrees):
    """ Converts polar coordinates to cartesian coordinates

    """
    x = radius * math.cos(math.radians(degrees))
    y = radius * math.sin(math.radians(degrees))
    return x, y


def draw_population_structure_graph(output_file_name,
                                    data,
                                    genotypes_ordering,
                                    struct_plot_settings):
    """ Draws a structure plot from the splicing data, broken up by genotype

    output_file_name is the file path where the plot will be written
    data is a pandas.DataFrame object representing the splicing data
    width is the width of the plot
    height is the height of the plot
    left_margin is the amount of white space on the left side of the plot
    right_margin is the amount of white space on the right side of the plot
    bottom_margin is the amount of white space below the plot
    top_margin is the amount of white space above the plot
    colors is a list of RGB objects which determine the colors of the bars in the structure plot
    axis_color is an RGB object representing the color of the axes, labels, and tick marks
    axis_thickness is the thickness of the axes and ticks
    tick_length is the length of the tick marks
    horiz_label_size is the font size of the labels along the horizontal (x) axis
    horiz_label_spacing represents the distance from the horizontal axis labels are drawn below the horizontal axis
    horiz_axis_title_size is the font size of the label for the horizontal axis
    horiz_axis_title_spacing is the distance from the horizontal axis that the title is drawn below the horizontal axis
    use_vertical_ticks is a boolean determining whether or not tick marks are drawn for the vertical (y) axis
    vertical_tick_spacing represents the scale for the tick marks on the vertical axis. Must be between 0 and 1
    vert_label_size is the font size of the vertical axis tick mark labels
    vert_label_spacing represents the distance from which the vertical axis labels are drawn to the left of the vertical axis
    include_key is a boolean determining whether or not a key should be drawn
    key_position is an array representing the cartesian coordinates of the key. [0,0] refers to the bottom left corner of the plot
    key_font_size is the font size of each entry in the key
    key_text_color is an RGB object which represents the color of the text in the key
    """

    # parse settings
    width=struct_plot_settings['plot_width']
    height=struct_plot_settings['plot_height']
    left_margin=struct_plot_settings['left_margin']
    right_margin=struct_plot_settings['right_margin']
    bottom_margin=struct_plot_settings['bottom_margin']
    top_margin=struct_plot_settings['top_margin']
    colors=struct_plot_settings['colors']
    axis_color=struct_plot_settings['axis_color']
    axis_thickness =struct_plot_settings['axis_thickness']
    tick_length=struct_plot_settings['tick_length']
    horiz_label_size=struct_plot_settings['horiz_label_size']
    horiz_label_spacing=struct_plot_settings['horiz_label_spacing']
    horiz_axis_title_size=struct_plot_settings['horiz_axis_title_size']
    horiz_axis_title_spacing=struct_plot_settings['horiz_axis_title_spacing']
    use_vertical_ticks=struct_plot_settings['use_vertical_ticks']
    vertical_tick_spacing=struct_plot_settings['vertical_tick_spacing']
    vert_label_size=struct_plot_settings['vert_label_size']
    vert_label_spacing=struct_plot_settings['vert_label_spacing']
    include_key=struct_plot_settings['include_key']
    key_position=struct_plot_settings['key_position']
    key_font_size=struct_plot_settings['key_font_size']
    key_text_color=struct_plot_settings['key_text_color']


    # sort the data frame so that homozygous reference is first, heterozygous is in the middle, etc
    indiv_and_genotypes = data.iloc[:,0]
    indivs_by_genotype = {}
    for indiv_id, genotype in indiv_and_genotypes.iteritems():
        if genotype not in indivs_by_genotype:
            indivs_by_genotype[genotype] = []
        indivs_by_genotype[genotype].append(indiv_id)

    new_index_order = []
    for genotype in genotypes_ordering:
        new_index_order.extend(indivs_by_genotype[genotype])
    sorted_data = data.reindex(index=new_index_order)

    # assign colors for each junction
    color_lookup = {}
    all_junctions = list(sorted_data.columns[1:])
    for i in range(len(all_junctions)):
        color_lookup[all_junctions[i]] = colors[i]

    width_per_individual = (width - (left_margin + right_margin)) / len(sorted_data.index)

    
    # create svg file
    f1 = open(output_file_name,'w+')

    # set up svg
    f1.write('<?xml version="1.0"?>\n')
    f1.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
    f1.write('<svg width="{0}in" height="{1}in" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {0} {1}">\n'.format(width,height))

    # flip to cartesian coordinates, (0,0) is at bottom left corner
    f1.write('<g transform="translate(0,{0}) scale(1,-1)">\n'.format(height))

    # draw the rectangles. also draw vertical separators and horizontal axis labels
    previous_width = left_margin
    current_width = left_margin
    current_genotype = sorted_data.iloc[0,0]
    for row_number in range(sorted_data.shape[0]):
        individual = sorted_data.iloc[row_number]

        # draw rectangles
        expression_total_value = 0
        for junction in all_junctions:

            junction_color = color_lookup[junction]

            current_junction_expression = individual[junction]
            start_y, rect_height = (height - top_margin - bottom_margin) * expression_total_value + bottom_margin, (height - top_margin - bottom_margin) * current_junction_expression            
            drawing.draw_rectangle(f1,current_width,start_y,width_per_individual,rect_height,junction_color)
            expression_total_value += current_junction_expression

        # draw separators and horizontal axis labels
        if individual.iloc[0] != current_genotype or row_number == sorted_data.shape[0] - 1:
            if individual.iloc[0] != current_genotype:
                drawing.draw_line(f1,current_width,bottom_margin-tick_length*0.5,current_width,height-top_margin,axis_thickness,axis_color)

            if horiz_label_size:
                drawing.draw_text(f1,current_genotype,previous_width+(current_width-previous_width)*0.5,bottom_margin-horiz_label_spacing,horiz_label_size,angle=0,color=axis_color)

            previous_width = current_width
            current_genotype = individual.iloc[0]

        current_width += width_per_individual


    # draw vertical tick marks
    if use_vertical_ticks:
        i = vertical_tick_spacing
        while i <= 1:
            y_coordinate = (height - top_margin - bottom_margin) * i + bottom_margin
            x_coordinate = left_margin - vert_label_spacing

            drawing.draw_line(f1,left_margin-tick_length*0.5,y_coordinate,left_margin+tick_length*0.5,y_coordinate,axis_thickness,axis_color)
            drawing.draw_text(f1,i,x_coordinate,y_coordinate-vert_label_size*0.25,vert_label_size,angle=0,color=axis_color)

            drawing.draw_line
            i += vertical_tick_spacing


    # draw the plot border
    left_x, right_x = left_margin, current_width
    top_y, bottom_y = height - top_margin, bottom_margin

    drawing.draw_line(f1,left_x-axis_thickness/2,top_y,right_x+axis_thickness/2,top_y,axis_thickness,axis_color)
    drawing.draw_line(f1,left_x,top_y+axis_thickness/2,left_x,bottom_y-axis_thickness/2,axis_thickness,axis_color)
    drawing.draw_line(f1,right_x,bottom_y-axis_thickness/2,right_x,top_y+axis_thickness/2,axis_thickness,axis_color)
    drawing.draw_line(f1,right_x+axis_thickness/2,bottom_y,left_x-axis_thickness/2,bottom_y,axis_thickness,axis_color)


    # draw horizontal axis title, if desired
    if horiz_axis_title_size != 0:
        drawing.draw_text(f1,data.columns[0],(left_margin+current_width)*0.5,bottom_margin-horiz_axis_title_spacing,horiz_axis_title_size,angle=0,color=axis_color)


    # draw the key if desired
    if include_key:

        spacing = key_font_size * 2

        for i in range(len(all_junctions)):
            junction_name = all_junctions[i]
            color = color_lookup[junction_name]

            drawing.draw_rectangle(f1,key_position[0],key_position[1] - i * spacing,key_font_size,key_font_size,color)
            drawing.draw_text_left(f1,data.columns[i+1],key_position[0]+spacing,key_position[1] - i * spacing,key_font_size,color=key_text_color)
    
            
    f1.write('</g>\n')
    f1.write('</svg>\n')

    f1.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw a hive plot from the given data')
    parser.add_argument('-s', '--settings', help='settings file, see sample1/test_settings for example', type=str, default = None)

    args = parser.parse_args()
    data, hive_plot_settings, struct_plot_settings = plot_settings.parse_settings(args.settings)
    

    if hive_plot_settings['draw_hive_plot']:
        print 'Drawing hive plot...'

        plotter_svg(file_name=hive_plot_settings['output_file_path'],
                    data=data,
                    hive_plot_settings=hive_plot_settings)

    if struct_plot_settings['draw_struct_plot']:

        print 'Drawing structure plot...'

        draw_population_structure_graph(output_file_name=struct_plot_settings['output_file_path'],
                                        data=data,
                                        struct_plot_settings=struct_plot_settings)

    print 'Done!'
