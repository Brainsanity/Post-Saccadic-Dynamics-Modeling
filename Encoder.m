classdef Encoder < handle

	properties (SetAccess = private)
		layers;%(4) = struct( 'name', [], 'locations', [], 'sRFRc', [], 'sRFKc', [], 'sRFRs', [], 'sRFKs', [] );
	end

	methods
		function obj = Encoder()
			obj.layers = struct( 'name', cell(1,4), 'locations', [], 'sRFRc', [], 'sRFKc', [], 'sRFRs', [], 'sRFKs', [] );
			names = {'POn', 'POff', 'MOn', 'MOff'};
			[obj.layers.name] = names{:};
			for( k =  1 : size(obj.layers,2) )
				obj.layers(k).locations = load( ['../../Data/Mosaic_' names{k} '_Radius20.0deg_maxMovPrctile20.mat'], 'rfPositions' );
				obj.layers(k).locations = obj.layers(k).locations.rfPositions;
			end
		end


		function DisplayMosaics(obj, radius)
			% figure( 'NumberTitle', 'off', 'name', 'Mosaics', 'color', 'k' );
			colors = {'r', 'b', 'm', 'c'};
			for( k = 1 : size(obj.layers,2) )
				figure( 'NumberTitle', 'off', 'name', sprintf( 'Mosaics | %s | [%.1f, %.1f]', obj.layers(k).name], -radius, radius ), 'color', 'k' );
				% subplot(2,2,k);
				plot( obj.layers(k).locations(:,1), obj.layers(k).locations(:,2), '.', 'color', colors{k}, 'markersize', 4 );
				% title( obj.layers(k).name );
				text( 0, radius*1.02, obj.layers(k).name, 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom', 'fontsize', 20, 'color', 'w' );
				axis equal;
				set( gca, 'xlim', [-radius radius], 'ylim', [-radius radius], 'fontsize', 20, 'visible', 'off', 'color', 'k', 'xcolor', 'w', 'ycolor', 'w' );
			end
		end
	end
end