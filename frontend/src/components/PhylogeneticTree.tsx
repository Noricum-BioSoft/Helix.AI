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
  
  if (!newick) {
    console.log('No newick data provided');
    return <div className="alert alert-warning">No tree data provided.</div>;
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
      const treeData = parseNewick(newick);
      console.log('Parsed tree data:', treeData);
      
      // Create hierarchy
      const hierarchy = d3.hierarchy(treeData);
      
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
      
      // Set up the tree layout with dynamic dimensions based on label length
      const baseWidth = 600;
      const labelWidth = Math.max(maxLabelLength * 8, 200); // 8px per character, minimum 200px
      const width = baseWidth + labelWidth;
      const height = 400;
      const treeLayout = d3.tree()
        .size([height, width - labelWidth - 50]); // Leave room for labels on the right
      
      const tree = treeLayout(hierarchy as any);
      
      // Update SVG width to accommodate labels
      svgRef.current!.setAttribute('width', String(width + 100));
      
      // Create SVG
      const svg = d3.select(svgRef.current);
      const g = svg.append('g')
        .attr('transform', 'translate(50, 50)');
      
      // Add links with error checking
      g.selectAll('.link')
        .data(tree.links())
        .enter()
        .append('path')
        .attr('class', 'link')
        .attr('fill', 'none')
        .attr('stroke', '#555')
        .attr('stroke-width', 1)
        .attr('d', (d: any) => {
          // Check for valid coordinates
          if (isNaN(d.source.x) || isNaN(d.source.y) || isNaN(d.target.x) || isNaN(d.target.y)) {
            console.error('Invalid coordinates in tree:', d);
            return '';
          }
          return d3.linkHorizontal()
            .x((d: any) => d.y)
            .y((d: any) => d.x)(d as any);
        });
      
      // Add nodes with error checking
      const node = g.selectAll('.node')
        .data(tree.descendants())
        .enter()
        .append('g')
        .attr('class', 'node')
        .attr('transform', (d: any) => {
          // Check for valid coordinates
          if (isNaN(d.x) || isNaN(d.y)) {
            console.error('Invalid node coordinates:', d);
            return 'translate(0,0)';
          }
          return `translate(${d.y},${d.x})`;
        });
      
      // Add circles for nodes
      node.append('circle')
        .attr('r', 3)
        .attr('fill', (d: any) => d.children ? '#555' : '#69b3a2');
      
      // Add labels - show full names, no truncation
      const labels = node.append('text')
        .attr('dy', '.31em')
        .attr('x', (d: any) => d.children ? -8 : 8)
        .attr('text-anchor', (d: any) => d.children ? 'end' : 'start')
        .text((d: any) => d.data.name || '')
        .style('font-size', '11px')
        .style('font-family', 'monospace')
        .style('cursor', 'pointer');
      
      // Add tooltips to show full names on hover (helpful for very long names)
      labels.append('title')
        .text((d: any) => d.data.name || '');
        
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
      <div style={{ height: '500px', width: '100%', overflow: 'auto' }}>
        <svg
          ref={svgRef}
          width="800"
          height="500"
          style={{ border: '1px solid #ccc', minWidth: '800px' }}
        />
      </div>
    </div>
  );
}; 