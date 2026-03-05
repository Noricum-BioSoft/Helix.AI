import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';

interface PhylogeneticTreeProps {
  newick: string;
}

interface TreeNode {
  name: string;
  children?: TreeNode[];
  branchLength?: number;
}

export const PhylogeneticTree: React.FC<PhylogeneticTreeProps> = ({ newick }) => {
  const svgRef = useRef<SVGSVGElement>(null);

  console.log('PhylogeneticTree component rendered with newick:', newick);
  console.log('Newick type:', typeof newick);
  console.log('Newick length:', newick?.length);
  console.log('Newick preview:', newick?.substring(0, 100));
  
  if (!newick) {
    console.log('No newick data provided');
    return <div className="alert alert-warning">No tree data provided.</div>;
  }
  
  if (typeof newick !== 'string') {
    console.error('Newick is not a string:', typeof newick, newick);
    return <div className="alert alert-danger">Invalid tree data format. Expected string, got {typeof newick}.</div>;
  }

  // Improved Newick parser
  const parseNewick = (newick: string): TreeNode => {
    try {
      // Remove semicolon at the end
      const cleanNewick = newick.replace(/;$/, '');
      
      const parseNode = (str: string): TreeNode => {
        str = str.trim();
        
        // If it's a leaf node (no parentheses)
        if (!str.includes('(')) {
          const match = str.match(/^([^:]+)(?::([^:]+))?$/);
          if (match) {
            return {
              name: match[1],
              branchLength: match[2] ? parseFloat(match[2]) : 0
            };
          }
          return { name: str };
        }
        
        // Find the matching parentheses
        let depth = 0;
        let start = str.indexOf('(');
        let end = -1;
        
        for (let i = start; i < str.length; i++) {
          if (str[i] === '(') depth++;
          else if (str[i] === ')') {
            depth--;
            if (depth === 0) {
              end = i;
              break;
            }
          }
        }
        
        if (end === -1) {
          console.error('Unmatched parentheses in Newick string');
          return { name: 'Error' };
        }
        
        // Extract the content inside parentheses
        const innerContent = str.substring(start + 1, end);
        const afterContent = str.substring(end + 1);
        
        // Parse children by splitting on commas at the top level
        const children: TreeNode[] = [];
        let currentChild = '';
        let childDepth = 0;
        
        for (let i = 0; i < innerContent.length; i++) {
          const char = innerContent[i];
          if (char === '(') childDepth++;
          else if (char === ')') childDepth--;
          else if (char === ',' && childDepth === 0) {
            if (currentChild.trim()) {
              children.push(parseNode(currentChild.trim()));
            }
            currentChild = '';
            continue;
          }
          currentChild += char;
        }
        
        if (currentChild.trim()) {
          children.push(parseNode(currentChild.trim()));
        }
        
        // Parse the name and branch length after the closing parenthesis
        const nameMatch = afterContent.match(/^([^:]*)(?::([^:]*))?$/);
        const name = nameMatch ? nameMatch[1] || 'Internal' : 'Internal';
        const branchLength = nameMatch && nameMatch[2] ? parseFloat(nameMatch[2]) : 0;
        
        return {
          name,
          children,
          branchLength
        };
      };
      
      return parseNode(cleanNewick);
    } catch (error) {
      console.error('Error parsing Newick string:', error);
      return { name: 'Parse Error' };
    }
  };

  useEffect(() => {
    if (!svgRef.current || !newick) return;

    try {
      // Clear previous content
      d3.select(svgRef.current).selectAll("*").remove();
      
      // Parse the Newick string
      console.log('🔍 Parsing Newick string, length:', newick.length);
      console.log('🔍 Newick string:', newick);
      const treeData = parseNewick(newick);
      console.log('🔍 Parsed tree data:', JSON.stringify(treeData, null, 2));
      
      // Validate parsed tree
      if (!treeData || !treeData.children || treeData.children.length === 0) {
        console.error('❌ Parsed tree has no children:', treeData);
        throw new Error('Failed to parse tree structure - no children found');
      }
      
      // Create hierarchy
      const hierarchy = d3.hierarchy(treeData);
      console.log('🔍 Created hierarchy, nodes:', hierarchy.descendants().length);
      
      // Validate hierarchy
      if (!hierarchy || hierarchy.descendants().length === 0) {
        throw new Error('Hierarchy has no nodes');
      }
      
      // Calculate max label length from hierarchy data
      const getAllNames = (node: any): string[] => {
        const names = [node.data.name || ''];
        if (node.children) {
          node.children.forEach((child: any) => {
            names.push(...getAllNames(child));
          });
        }
        return names;
      };
      const allNames = getAllNames(hierarchy);
      const maxLabelLength = Math.max(...allNames.map(name => name.length), 10);
      console.log('🔍 Max label length:', maxLabelLength);
      
      // Set up the tree layout with dynamic dimensions based on label length
      const baseWidth = 600;
      const labelWidth = Math.max(maxLabelLength * 8, 200); // 8px per character, minimum 200px
      const width = baseWidth + labelWidth;
      const height = Math.max(400, hierarchy.leaves().length * 50); // Dynamic height based on number of leaves
      console.log('🔍 Tree dimensions:', { width, height, labelWidth });
      
      const treeLayout = d3.tree()
        .size([height - 100, width - labelWidth - 50]); // Leave room for labels on the right
      
      const tree = treeLayout(hierarchy as any);
      console.log('🔍 Tree layout computed, links:', tree.links().length, 'nodes:', tree.descendants().length);
      
      // Update SVG width to accommodate labels
      const svgWidth = width + 100;
      const svgHeight = height + 100;
      svgRef.current!.setAttribute('width', String(svgWidth));
      svgRef.current!.setAttribute('height', String(svgHeight));
      console.log('🔍 SVG dimensions set:', { svgWidth, svgHeight });
      
      // Create SVG - ensure it's cleared first
      const svg = d3.select(svgRef.current);
      if (!svgRef.current) {
        throw new Error('SVG ref is null');
      }
      
      // Clear any existing content
      svg.selectAll("*").remove();
      console.log('🔍 SVG cleared');
      
      // Set SVG background
      svg.attr('style', 'background-color: #fff; display: block;');
      
      const g = svg.append('g')
        .attr('transform', 'translate(50, 50)');
      console.log('🔍 SVG group created');
      
      // Add links with error checking
      const links = tree.links();
      console.log('🔍 Adding', links.length, 'links');
      const linkSelection = g.selectAll('.link')
        .data(links)
        .enter()
        .append('path')
        .attr('class', 'link')
        .attr('fill', 'none')
        .attr('stroke', '#555')
        .attr('stroke-width', 2);
      
      linkSelection.attr('d', (d: any) => {
        // Check for valid coordinates
        if (isNaN(d.source.x) || isNaN(d.source.y) || isNaN(d.target.x) || isNaN(d.target.y)) {
          console.error('❌ Invalid coordinates in tree link:', d);
          return '';
        }
        const linkGenerator = d3.linkHorizontal()
          .x((d: any) => d.y)
          .y((d: any) => d.x);
        const path = linkGenerator(d as any);
        console.log('🔍 Link path:', path, 'from', d.source.data.name, 'to', d.target.data.name);
        return path;
      });
      console.log('🔍 Links added');
      
      // Add nodes with error checking
      const nodes = tree.descendants();
      console.log('🔍 Adding', nodes.length, 'nodes');
      const nodeSelection = g.selectAll('.node')
        .data(nodes)
        .enter()
        .append('g')
        .attr('class', 'node')
        .attr('transform', (d: any) => {
          // Check for valid coordinates
          if (isNaN(d.x) || isNaN(d.y)) {
            console.error('❌ Invalid node coordinates:', d);
            return 'translate(0,0)';
          }
          const transform = `translate(${d.y},${d.x})`;
          console.log('🔍 Node transform:', transform, 'for', d.data.name);
          return transform;
        });
      
      // Add circles for nodes
      nodeSelection.append('circle')
        .attr('r', 4)
        .attr('fill', (d: any) => {
          const color = d.children ? '#555' : '#69b3a2';
          console.log('🔍 Node color:', color, 'for', d.data.name, 'has children:', !!d.children);
          return color;
        })
        .attr('stroke', '#fff')
        .attr('stroke-width', 1);
      
      // Add labels - show full names, no truncation
      const labels = nodeSelection.append('text')
        .attr('dy', '.31em')
        .attr('x', (d: any) => {
          const offset = d.children ? -10 : 10;
          console.log('🔍 Label offset:', offset, 'for', d.data.name);
          return offset;
        })
        .attr('text-anchor', (d: any) => d.children ? 'end' : 'start')
        .text((d: any) => {
          const name = d.data.name || '';
          console.log('🔍 Label text:', name);
          return name;
        })
        .style('font-size', '12px')
        .style('font-family', 'monospace')
        .style('fill', '#333')
        .style('cursor', 'pointer');
      
      // Add tooltips to show full names on hover (helpful for very long names)
      labels.append('title')
        .text((d: any) => d.data.name || '');
      
      console.log('✅ Tree rendering complete');
        
          } catch (error) {
        console.error('Error rendering phylogenetic tree:', error);
        // Show error message in the SVG
        d3.select(svgRef.current)
          .append('text')
          .attr('x', 50)
          .attr('y', 50)
          .attr('fill', 'red')
          .text('Error rendering tree: ' + (error instanceof Error ? error.message : String(error)));
      }
  }, [newick]);

  return (
    <div className="bg-light p-3 border rounded mb-3">
      <h5>Phylogenetic Tree</h5>
      {!newick && (
        <div className="alert alert-warning">No tree data available</div>
      )}
      {newick && (
        <div style={{ height: '600px', width: '100%', overflow: 'auto', backgroundColor: '#fff' }}>
          <svg
            ref={svgRef}
            width="800"
            height="600"
            style={{ 
              border: '1px solid #ccc', 
              minWidth: '800px',
              display: 'block',
              backgroundColor: '#fff'
            }}
          />
        </div>
      )}
    </div>
  );
}; 